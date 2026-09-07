// SNP-site detection straight from the imported hifiasm CIGARs.
//
// The CIGARs already record, for every overlap, exactly which columns disagree
// (hifiasm emits '=' and 'X' as distinct ops, preserved through the ingest
// transpose as CigarOpMatch / CigarOpMismatch). So finding candidate sites needs
// no alignment and no MSA -- only a walk over token runs, which is far cheaper
// than a per-base scan because match runs collapse to a single token. Measured
// on the 251k-overlap dataset: 112,979,378 tokens over ~2.8G aligned bases.
//
// The shape follows hifiasm's own error-correction design (Correct.cpp,
// markSNP / addSNPtohaplotype -- present but commented out in the vendored
// fork), which is two-pass for a good reason:
//
//   Pass 1  For each read r, walk every overlap involving r and count, per
//           position of r, how many partners disagree and how many cover it.
//           Cheap: two counters per position, no per-site bookkeeping.
//
//   Scan    Positions where enough partners disagree are candidate sites.
//
//   Pass 2  (not done here) Re-walk to collect which partner carries which
//           allele -- needed to build anchors, because the AGREEING reads are
//           as much a part of a site as the disagreeing ones.
//
// This file implements pass 1 and the scan, and reports the distribution. It
// creates no anchors and changes nothing downstream: the point is to find out
// whether the disagreement-fraction signal actually separates het sites from
// sequencing error before anything is built on top of it. The expected shape is
// bimodal -- a read's own error makes ~all partners disagree, a het site makes
// ~half disagree -- but ONT errors are systematic in homopolymers, so whether
// that separation is real is an empirical question, not an assumption.

#include "Assembler.hpp"
#include "HifiasmImportedCigarStore.hpp"
#include "OverlapCigarStore.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



void Assembler::detectCigarSnpSites(
    uint64_t minDisagreeCount,
    double minDisagreeFraction,
    double maxDisagreeFraction,
    uint64_t threadCount)
{
    if(hifiasmImportedCigarStore.empty()) {
        cout << timestamp << "SNP-site detection: the imported CIGAR store is "
            "empty (Align.useHifiasmBaseAlignment false?). Nothing to do." << endl;
        return;
    }
    const auto tBegin = steady_clock::now();
    if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
    if(threadCount == 0) threadCount = 1;

    const ReadId readCount = ReadId(reads->readCount());

    // Index the store by read. It is keyed by (pairKey, strand) for per-pair
    // lookup, so "every overlap involving read r" needs building once here.
    // Records whose CIGAR does not span its declared box are dropped: hifiasm
    // concatenates only the windows that aligned, so a dropped window leaves the
    // tokens describing a shorter span than the coordinates claim, and the token
    // stream carries no marker of where the hole was. Walking such a record from
    // qStart/tStart silently desynchronises (measured: 0.12% of records on the
    // GIAB fixture, 0.48% at 251k overlaps).
    vector<vector<const HifiasmImportedCigarStore::Record*>> byRead(readCount);
    uint64_t usable = 0, skippedShort = 0, skippedNoCigar = 0;
    hifiasmImportedCigarStore.forEachRecord(
        [&](const HifiasmImportedCigarStore::Record& rec) {
            if(rec.cigarTokenCount == 0) { ++skippedNoCigar; return; }
            if(rec.cigarQuerySpan  != rec.qEnd - rec.qStart ||
               rec.cigarTargetSpan != rec.tEnd - rec.tStart) { ++skippedShort; return; }
            if(rec.readIdQ >= readCount || rec.readIdT >= readCount) return;
            byRead[rec.readIdQ].push_back(&rec);
            byRead[rec.readIdT].push_back(&rec);
            ++usable;
        });

    cout << timestamp << "SNP-site detection: " << usable << " usable records ("
         << skippedShort << " skipped, CIGAR shorter than its box; "
         << skippedNoCigar << " with no CIGAR)." << endl;

    // Disagreement-fraction histogram, in twentieths, plus the site tally.
    std::array<std::atomic<uint64_t>, 21> fractionHistogram{};
    for(auto& bucket: fractionHistogram) bucket.store(0);
    std::atomic<uint64_t> candidateSites{0}, positionsExamined{0};

    std::atomic<ReadId> nextReadId{0};
    auto worker = [&]() {
        vector<uint16_t> disagree;
        vector<int32_t> coverDelta;
        for(;;) {
            const ReadId readId = nextReadId.fetch_add(1);
            if(readId >= readCount) break;
            const auto& records = byRead[readId];
            if(records.empty()) continue;

            const uint32_t readLength = uint32_t(reads->getRead(readId).baseCount);
            disagree.assign(readLength, 0);
            // Coverage as a DIFFERENCE array: each record covers one contiguous
            // interval of this read, so it costs two increments per record and a
            // single prefix sum at the end -- not a per-base loop. This is what
            // lets the token walk stay O(mismatches) rather than O(aligned
            // bases); see the match-run skip below.
            coverDelta.assign(size_t(readLength) + 1, 0);

            for(const HifiasmImportedCigarStore::Record* rec: records) {
                const bool selfIsQuery = (rec->readIdQ == readId);
                // The tokens run query-forward and target-in-alignment-
                // orientation. When this read is the target of a reverse-strand
                // overlap its alignment-orientation position a maps to forward
                // position (len - 1 - a); tLen - tEnd is where that frame starts.
                const uint32_t targetLength =
                    uint32_t(reads->getRead(ReadId(rec->readIdT)).baseCount);
                uint64_t qPos = rec->qStart;
                uint64_t tPos = rec->isSameStrand ?
                    rec->tStart : (targetLength - rec->tEnd);

                // Coverage of this record on this read: one interval, recorded
                // in the difference array. Both endpoints are already known
                // from the box, so no walking is needed for it.
                {
                    const uint32_t from = selfIsQuery ? rec->qStart : rec->tStart;
                    const uint32_t to   = selfIsQuery ? rec->qEnd   : rec->tEnd;
                    if(to > from && to <= readLength) {
                        coverDelta[from]++;
                        coverDelta[to]--;
                    }
                }

                // Now walk only for the disagreements. Following hifiasm's
                // markSNP_detail: a match run advances by its whole length in
                // one step -- matches are the overwhelming majority of columns,
                // so looping over them would dominate everything. Only mismatch
                // runs are visited per base. Indel runs consume one side only
                // and say nothing about a position's allele, so they just
                // advance the relevant cursor.
                for(const CigarToken token: hifiasmImportedCigarStore.tokensOf(*rec)) {
                    const uint8_t op = token.op();
                    const uint16_t length = token.len();
                    if(op == CigarOpMismatch) {
                        for(uint16_t i = 0; i < length; i++) {
                            uint64_t self;
                            if(selfIsQuery) {
                                self = qPos + i;
                            } else if(rec->isSameStrand) {
                                self = tPos + i;
                            } else {
                                self = uint64_t(targetLength) - 1 - (tPos + i);
                            }
                            if(self < readLength && disagree[self] < 0xffff) {
                                disagree[self]++;
                            }
                        }
                    }
                    if(opConsumesQuery(op))  qPos += length;
                    if(opConsumesTarget(op)) tPos += length;
                }
            }

            // Scan. A position is a candidate when enough partners disagree AND
            // the disagreeing share looks like a second haplotype rather than
            // this read's own error: near 1.0 means every partner disagrees, so
            // the odd base is almost certainly ours.
            uint64_t localSites = 0, localExamined = 0;
            std::array<uint64_t, 21> localHistogram{};
            int32_t running = 0;
            for(uint32_t position = 0; position < readLength; position++) {
                running += coverDelta[position];      // prefix sum -> coverage
                const uint32_t total = uint32_t(running < 0 ? 0 : running);
                if(total == 0) continue;
                ++localExamined;
                const uint32_t bad = disagree[position];
                if(bad == 0) continue;
                const double fraction = double(bad) / double(total);
                localHistogram[std::min<size_t>(20, size_t(fraction * 20.0))]++;
                if(bad >= minDisagreeCount &&
                   fraction >= minDisagreeFraction && fraction <= maxDisagreeFraction) {
                    ++localSites;
                }
            }
            candidateSites.fetch_add(localSites, std::memory_order_relaxed);
            positionsExamined.fetch_add(localExamined, std::memory_order_relaxed);
            for(size_t i = 0; i < localHistogram.size(); i++) {
                if(localHistogram[i]) {
                    fractionHistogram[i].fetch_add(localHistogram[i],
                        std::memory_order_relaxed);
                }
            }
        }
    };

    vector<std::thread> threads;
    for(uint64_t i = 0; i < threadCount; i++) threads.emplace_back(worker);
    for(auto& t: threads) t.join();

    const double elapsed = seconds(steady_clock::now() - tBegin);
    uint64_t withAnyDisagreement = 0;
    for(const auto& bucket: fractionHistogram) withAnyDisagreement += bucket.load();

    cout << timestamp << "SNP-site detection: examined "
         << positionsExamined.load() << " covered positions, "
         << withAnyDisagreement << " with at least one disagreeing partner, "
         << candidateSites.load() << " candidate sites (>= " << minDisagreeCount
         << " disagreeing, fraction in [" << minDisagreeFraction << ", "
         << maxDisagreeFraction << "]), in " << elapsed << " s." << endl;
    cout << "  disagreement-fraction histogram (bucket = fraction of covering "
            "partners that disagree):" << endl;
    for(size_t i = 0; i < fractionHistogram.size(); i++) {
        const uint64_t count = fractionHistogram[i].load();
        if(count == 0) continue;
        cout << "    [" << (0.05 * double(i)) << " - " << (0.05 * double(i + 1))
             << "): " << count << endl;
    }
    performanceLog << timestamp << "SNP-site detection took " << elapsed << " s." << endl;
}
