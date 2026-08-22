#include "Shasta2SnpmerAnchors.hpp"

#include "Assembler.hpp"
#include "Shasta2Anchors.hpp"
#include "Reads.hpp"
#include "ReadId.hpp"
#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "timestamp.hpp"
#include "dinaraTypes.hpp"
#include "span.hpp"
#include "DINARA_ASSERT.hpp"

// Bundled myloasm SNPmer engine C ABI (via the hifiasm_overlaps static lib).
#include "myloasm_ffi.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <unordered_map>

using namespace dinara;
using std::cout;
using std::endl;
using std::string;
using std::vector;



uint64_t dinara::appendSnpmerHetAnchors(
    Assembler& assembler,
    const vector<string>& inputFileNames,
    uint64_t threadCount,
    uint64_t minMembersPerAllele,
    bool requireBiallelic)
{
    const Reads& reads = assembler.getReads();
    const auto anchors = assembler.shasta2Anchors;   // shared_ptr<Shasta2Anchors>
    DINARA_ASSERT(anchors);

    // ------------------------------------------------------------------ //
    // 1. Index the input reads with myloasm (SNPmer + syncmer markers).   //
    //    myloasm re-parses the files; it may index a different subset than //
    //    dinara's store, so members are reconciled BY NAME below.          //
    // ------------------------------------------------------------------ //
    vector<const char*> paths;
    paths.reserve(inputFileNames.size());
    for(const string& f : inputFileNames) {
        paths.push_back(f.c_str());
    }

    MyloReadIndex idx;
    std::memset(&idx, 0, sizeof(idx));
    const int rc = myloasm_index_reads(
        paths.data(), paths.size(),
        /*kmer_size*/ 0, /*c*/ 0, int(threadCount), &idx);
    if(rc != 0) {
        cout << timestamp << "  myloasm_index_reads failed (rc=" << rc
             << "); no snpmer anchors appended." << endl;
        return 0;
    }
    const int k = idx.k;
    const int mid = (k - 1) / 2;
    cout << timestamp << "  myloasm indexed " << idx.n_reads
         << " reads, k=" << k << "." << endl;

    // ------------------------------------------------------------------ //
    // 2. Bucket every SNPmer occurrence by its CANONICAL key (= one allele //
    //    of one site) into (OrientedReadId, rawPosition) members. Track the //
    //    split (key with the middle base cleared) so the two alleles of a   //
    //    site group.                                                        //
    //                                                                       //
    // Every member is placed in the read's OWN sequencing frame (strand 0)  //
    // at rawPosition = snpCol - 1; appendHetAnchorPair generates the        //
    // reverse-complement (strand 1) placement itself. We deliberately do    //
    // NOT pre-orient members by a per-SNP canonical strand: a read has one  //
    // physical orientation, but two SNPs on it can each have their          //
    // canonical allele readable on OPPOSITE strands. Assigning per-SNP      //
    // strands puts both a strand-0 and a strand-1 member on one read; the   //
    // RC of the strand-1 member then lands on strand 0 on top of a          //
    // neighboring strand-0 canonical member, producing two anchors at the   //
    // same oriented-read position and crashing shasta2's journey builder.   //
    // Allele separation does not need orientation -- it comes entirely from //
    // the strand-invariant canonical key (allele0 and allele1 have distinct //
    // keys). The het-bubble engine likewise uses each read's single         //
    // alignment frame, never a per-SNP one.                                 //
    // ------------------------------------------------------------------ //
    struct Member { uint32_t orientedValue; uint32_t rawPosition; };
    struct Allele { vector<Member> members; };
    // key -> allele bucket
    std::unordered_map<uint64_t, Allele> alleleByKey;
    // split -> the (up to two) keys seen for it
    std::unordered_map<uint64_t, vector<uint64_t>> keysBySplit;

    // Middle base occupies bits [2*mid, 2*mid+1] of the packed k-mer.
    const uint64_t fullMask  = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    const uint64_t splitMask = fullMask & ~(3ULL << (2 * mid));

    uint64_t occTotal = 0, occUnmappedRead = 0, occOob = 0, occAmbig = 0,
             occKept = 0;

    for(size_t r = 0; r < idx.n_reads; ++r) {
        const MyloReadMarkers& R = idx.reads[r];

        // Reconcile myloasm read -> dinara ReadId by name (first token).
        const span<const char> nameSpan(R.name, R.name + R.name_len);
        const ReadId readId = reads.getReadId(nameSpan);
        if(readId == invalidReadId) {
            occUnmappedRead += R.n_snpmers;
            continue;
        }
        const uint32_t L = uint32_t(reads.getReadRawSequenceLength(readId));

        for(size_t s = 0; s < R.n_snpmers; ++s) {
            ++occTotal;
            const uint32_t pos = R.snpmers[s].pos;
            const uint64_t key = R.snpmers[s].key;

            // Need the full k-mer window around the SNP column.
            if(uint64_t(pos) + uint64_t(k) > L) { ++occOob; continue; }

            const uint32_t snpCol = pos + uint32_t(mid);   // forward SNP column

            // Place the k=2 het marker in the read's own (strand 0) frame:
            // predBase = snpCol-1, alleleBase = snpCol. appendHetAnchorPair
            // stores it at rawPosition + hetK/2 = snpCol and mirrors it onto
            // strand 1 itself. No per-SNP orientation (see the block comment
            // above for why that is unsafe).
            if(snpCol < 1) { ++occOob; continue; }
            const uint32_t rawPosition = snpCol - 1;
            if(uint64_t(rawPosition) + 1 >= L) { ++occOob; continue; }

            const OrientedReadId oid(readId, Strand(0));
            alleleByKey[key].members.push_back(
                Member{ oid.getValue(), rawPosition });

            const uint64_t split = key & splitMask;
            auto& keys = keysBySplit[split];
            if(std::find(keys.begin(), keys.end(), key) == keys.end()) {
                keys.push_back(key);
            }
            ++occKept;
        }
    }
    (void)occAmbig;

    // ------------------------------------------------------------------ //
    // 3. Per allele: dedup a read that occurs twice in the SAME allele     //
    //    (keep one member). Per split: drop any read appearing in BOTH      //
    //    alleles (paralog/repeat -- it cannot phase the site).              //
    // ------------------------------------------------------------------ //
    // Dedup within each allele by ReadId (both strands map to one ReadId).
    for(auto& kv : alleleByKey) {
        auto& m = kv.second.members;
        std::sort(m.begin(), m.end(),
            [](const Member& a, const Member& b) {
                return a.orientedValue < b.orientedValue;
            });
        m.erase(std::unique(m.begin(), m.end(),
            [](const Member& a, const Member& b) {
                return (a.orientedValue >> 1) == (b.orientedValue >> 1);
            }), m.end());
    }

    uint64_t paralogReadsDropped = 0;
    for(auto& sk : keysBySplit) {
        const vector<uint64_t>& keys = sk.second;
        if(keys.size() < 2) continue;
        // Count ReadId occurrences across the split's alleles.
        std::unordered_map<uint32_t, uint32_t> readCount;
        for(uint64_t key : keys) {
            for(const Member& m : alleleByKey[key].members) {
                readCount[m.orientedValue >> 1]++;
            }
        }
        // Remove reads that appear in more than one allele.
        for(uint64_t key : keys) {
            auto& m = alleleByKey[key].members;
            const size_t before = m.size();
            m.erase(std::remove_if(m.begin(), m.end(),
                [&](const Member& x) {
                    return readCount[x.orientedValue >> 1] > 1;
                }), m.end());
            paralogReadsDropped += (before - m.size());
        }
    }

    // ------------------------------------------------------------------ //
    // 4. Append surviving alleles as k=2 het anchors. Sites are emitted    //
    //    only when they still form a real bubble (>=1 allele, or both when  //
    //    requireBiallelic) and each kept allele has >= minMembersPerAllele. //
    // ------------------------------------------------------------------ //
    uint64_t sitesEmitted = 0, allelesAppended = 0, allelesSkippedSmall = 0,
             sitesSkippedNotBiallelic = 0;

    for(auto& sk : keysBySplit) {
        const vector<uint64_t>& keys = sk.second;

        // How many alleles of this split survive the member-count gate?
        uint64_t survivingAlleles = 0;
        for(uint64_t key : keys) {
            if(alleleByKey[key].members.size() >= minMembersPerAllele) {
                ++survivingAlleles;
            }
        }
        if(requireBiallelic && survivingAlleles < 2) {
            ++sitesSkippedNotBiallelic;
            continue;
        }
        if(survivingAlleles == 0) continue;

        bool emittedHere = false;
        for(uint64_t key : keys) {
            auto& m = alleleByKey[key].members;
            if(m.size() < minMembersPerAllele) {
                if(m.size() > 0) ++allelesSkippedSmall;
                continue;
            }
            vector<std::pair<OrientedReadId, uint32_t>> members;
            members.reserve(m.size());
            for(const Member& x : m) {
                members.push_back({ OrientedReadId::fromValue(x.orientedValue),
                                    x.rawPosition });
            }
            anchors->appendHetAnchorPair(members);
            ++allelesAppended;
            emittedHere = true;
        }
        if(emittedHere) ++sitesEmitted;
    }

    cout << timestamp << "  snpmer occurrences: " << occTotal
         << " total, " << occKept << " placed, "
         << occUnmappedRead << " on unmapped reads, "
         << occOob << " out-of-bounds." << endl;
    cout << timestamp << "  snpmer filtering: " << paralogReadsDropped
         << " paralog members dropped, " << allelesSkippedSmall
         << " alleles below minMembers(" << minMembersPerAllele << "), "
         << sitesSkippedNotBiallelic << " sites not biallelic." << endl;
    cout << timestamp << "  snpmer anchors appended: " << allelesAppended
         << " allele anchors across " << sitesEmitted << " sites." << endl;

    myloasm_read_index_free(&idx);
    return allelesAppended;
}
