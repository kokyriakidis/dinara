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
//   Pass 2  Re-walk, for OWNED sites only, to collect which partner carries
//           which allele -- needed to build anchors, because the AGREEING reads
//           are as much a part of a site as the disagreeing ones. "Owned" means
//           this read has the lowest ReadId among those covering the site, which
//           deduplicates the ~coverage-fold redundant detections before the
//           expensive pass rather than after it.
//
// MEASURED ACCURACY (989-read chr1:15.0-15.4Mb fixture, against the published
// marbl/HG002 mat-vs-pat het track; harness: evaluate2.py alongside the truth
// files, driven by DINARA_SNP_OWNED_DUMP / DINARA_SNP_SITE_DUMP):
//
//     precision 98.0%   recall 96.1%
//
// The DENOMINATOR matters more than the numbers, and getting it wrong made this
// detector look far worse than it is. Of 303 truth SNVs in the read-covered
// region:
//   - 79 sit in a homopolymer in one frame or the other. The truth track is
//     itself a minimap2 alignment, and inside a homopolymer the placement of a
//     difference is ambiguous (the same change spells as a substitution or as a
//     shifted indel), so those calls can neither confirm nor refute us. They
//     belong in neither numerator nor denominator. 224 remain.
//   - 44 of those 224 are MONOALLELIC in the read pileup -- verified against
//     minimap2 directly, independent of dinara's overlap set. They cluster in
//     two haplotype-dropout regions of this read subset (60 maternal reads to 1
//     paternal in one, 54 paternal to 10 maternal in the other). No detector can
//     call a het site where one haplotype is absent from the reads. 180 remain.
//
// Against those 180 reachable variants: 173 found, 7 missed (3 to the cc depth
// floor, 2 to the STR gate, 2 never became candidates). Reported against the
// naive 303 the same run scores 58.7%, which is why recall figures here always
// state their denominator.
//
// The threshold surface is FLAT: sweeping the minor-allele fraction over
// 0.15-0.30 and cc over 6-18, with the STR gate on or off, moves F1 only
// between 85.5 and 87.2. A minor-allele fraction of 0.25 is optimal at every cc.
// Thresholds are not the lever here; do not spend effort retuning them.
//
// Detection is deliberately PERMISSIVE, matching hifiasm's own criterion: a
// position is a candidate when at least 2 covering partners disagree, with no
// frequency cutoff (Correct.cpp: snp_threshold = 1, tested as
// flag[i] > snp_threshold, and nothing more). Deciding whether a minority
// allele is real belongs downstream -- to a binomial test against the assumed
// error rate, and to homopolymer/repeat/strand-bias filters -- not to a
// frequency cutoff here.
//
// That was measured rather than assumed. A [0.2, 0.8] fraction window cut 3662
// sites to 518 on the GIAB fixture: it does remove real junk (2129 monoallelic
// sites down to 49), but it also discards 1027 of 1493 BIALLELIC sites, 69% of
// them. Those are the low-VAF tails -- a second allele carried by a few of ~34
// reads, or the owning read itself carrying the rare allele -- which is exactly
// the population a statistical test exists to adjudicate and a frequency cutoff
// cannot.
//
// The reported disagreement-fraction histogram is trimodal on this data: a mass
// near 0 (one partner's error), a bump at ~0.5 (a second haplotype), and a spike
// at 1.0 (this read's own error). It is kept as a diagnostic, not used as a
// filter.
//
// This file creates no anchors and changes nothing downstream.

#include "Assembler.hpp"
#include "HifiasmImportedCigarStore.hpp"
#include "OverlapCigarStore.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"
#include "hetSignificance.hpp"
#include "PhasingKmeansTypes.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>
#include <vector>
#include <fstream>
#include <memory>
#include <sstream>
#include <cstdlib>

using namespace dinara;
using namespace std;



void Assembler::detectCigarSnpSites(
    uint64_t minDisagreeCount,
    double minDisagreeFraction,
    double maxDisagreeFraction,
    double hetErrorRate,
    double strandBiasPValue,
    double siteMinPurity,
    double siteMinAltDominance,
    double minAlleleFraction,
    bool filterHomopolymer,
    bool filterStr,
    double alleleCoverageRate,
    uint64_t alleleCoverageFloor,
    uint64_t ploidy,
    vector<CigarSnpSite>* sitesOut,
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

    // Dynamic minimum support for a het allele, derived from the dataset's own
    // coverage statistics rather than fixed. This is hifiasm's `cc`, from the
    // LIVE phasing path (gen_rphase_dp0_single_path, Correct.cpp:9435):
    //
    //     cc = ((het_cov > 0) ? het_cov : (hom_cov / n_hap));
    //     cc *= cut_rate; if(cc < cut_bd) cc = cut_bd;
    //
    // with the caller passing asm_opt.het_cov, asm_opt.hom_cov,
    // asm_opt.polyploidy, cut_rate = 0.7, cut_bd = 6. hifiasm gates on it
    // twice, both times against occ_0 -- the MAJOR allele's count -- so it
    // means "this site's dominant allele carries the read support a real
    // haplotype would". At a true het site each haplotype gets about het_cov
    // reads, so requiring 70% of that adapts to the data instead of guessing.
    //
    // dinara already has the same statistics: hifiasm's own filter reports its
    // coverage peaks and they are kept in kmerDistributionInfo for the run.
    const uint64_t coverageHet = assemblerInfo->kmerDistributionInfo.coverageHet;
    const uint64_t coverageHom = assemblerInfo->kmerDistributionInfo.coverageHom;
    uint64_t minAlleleCoverage = alleleCoverageFloor;
    {
        const uint64_t base = (coverageHet != invalid<uint64_t> && coverageHet > 0) ?
            coverageHet :
            ((coverageHom != invalid<uint64_t> && ploidy > 0) ? (coverageHom / ploidy) : 0);
        const uint64_t scaled = uint64_t(double(base) * alleleCoverageRate);
        if(scaled > minAlleleCoverage) minAlleleCoverage = scaled;
        cout << timestamp << "SNP-site detection: minimum major-allele support "
             << minAlleleCoverage << " (het coverage " << base << " x "
             << alleleCoverageRate << ", floored at " << alleleCoverageFloor
             << ")." << endl;
    }

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

    // DINARA_DUMP_IMPORTED_PAF: write the records exactly as imported, as a
    // PAF with CIGARs. This pins the comparison input: the reads, the filter
    // and the seeding are whatever THIS run used, so an independent
    // implementation reading this file sees precisely what detection sees.
    // Comparing against a separately-run hifiasm cannot do that -- its filter
    // and read set differ, and it legitimately spells alignments differently.
    //
    // Ops are written in DINARA's convention: '=' match, 'X' mismatch,
    // 'I' consumes the query, 'D' consumes the target. Column 1 is dinara's
    // query (hifiasm's x), NOT the PAF role order hifiasm's own writer uses.
    if(const char* pafPath = std::getenv("DINARA_DUMP_IMPORTED_PAF")) {
        ofstream paf(pafPath);
        static const char opChar[4] = {'=', 'X', 'I', 'D'};
        uint64_t written = 0;
        hifiasmImportedCigarStore.forEachRecord(
            [&](const HifiasmImportedCigarStore::Record& rec) {
                if(rec.cigarTokenCount == 0) return;
                if(rec.readIdQ >= readCount || rec.readIdT >= readCount) return;
                const auto qn = reads->getReadName(ReadId(rec.readIdQ));
                const auto tn = reads->getReadName(ReadId(rec.readIdT));
                paf.write(&*qn.begin(), std::streamsize(qn.size()));
                paf << '\t' << reads->getRead(ReadId(rec.readIdQ)).baseCount
                    << '\t' << rec.qStart << '\t' << rec.qEnd
                    << '\t' << (rec.isSameStrand ? '+' : '-') << '\t';
                paf.write(&*tn.begin(), std::streamsize(tn.size()));
                paf << '\t' << reads->getRead(ReadId(rec.readIdT)).baseCount
                    << '\t' << rec.tStart << '\t' << rec.tEnd
                    << "\t0\t0\t255\tcg:Z:";
                for(const CigarToken tk: hifiasmImportedCigarStore.tokensOf(rec)) {
                    paf << tk.len() << opChar[tk.op() & 3];
                }
                paf << '\n';
                ++written;
            });
        cout << timestamp << "SNP-site detection: wrote " << written
             << " imported records to " << pafPath << endl;
    }

    cout << timestamp << "SNP-site detection: " << usable << " usable records ("
         << skippedShort << " skipped, CIGAR shorter than its box; "
         << skippedNoCigar << " with no CIGAR)." << endl;

    // Disagreement-fraction histogram, in twentieths, plus the site tally.
    std::array<std::atomic<uint64_t>, 21> fractionHistogram{};
    for(auto& bucket: fractionHistogram) bucket.store(0);
    std::atomic<uint64_t> candidateSites{0}, positionsExamined{0};
    std::atomic<uint64_t> ownedSites{0}, allelesTotal{0}, significantSites{0};
    std::atomic<uint64_t> droppedHomopolymer{0}, droppedRepeat{0};
    std::atomic<uint64_t> droppedStrandBias{0}, sitesAfterFilters{0};
    std::atomic<uint64_t> droppedTiedAlleles{0}, droppedImpure{0}, droppedNotDominant{0};
    std::atomic<uint64_t> droppedLowCoverage{0}, droppedAdjacent{0};
    // Indexed [binomialPass*2 + ccPass], so cell 1 = passed cc but not the
    // binomial (= rejected by the binomial alone) and cell 2 = the reverse.
    // Cross-check when reading these: cells 0+2 must sum to droppedLowCoverage
    // (every cc rejection) and cells 0+1 to the binomial's own rejection count.
    std::array<std::atomic<uint64_t>, 4> binomialVsCc{};
    std::atomic<uint64_t> noAlternateSites{0};
    for(auto& cell: binomialVsCc) cell.store(0);
    std::mutex sitesMutex;
    std::array<std::atomic<uint64_t>, 5> significanceHistogram{};
    for(auto& bucket: significanceHistogram) bucket.store(0);
    std::atomic<uint64_t> matchColumnsChecked{0}, matchColumnsAgree{0};
    std::atomic<uint64_t> mismatchColumnsChecked{0}, mismatchColumnsDiffer{0};
    std::array<std::atomic<uint64_t>, 9> alleleCountHistogram{};
    for(auto& bucket: alleleCountHistogram) bucket.store(0);

    // Optional parity dump.
    std::unique_ptr<ofstream> dumpStream;
    std::mutex dumpMutex;
    if(const char* path = std::getenv("DINARA_SNP_DUMP")) {
        dumpStream = std::make_unique<ofstream>(path);
        cout << timestamp << "SNP-site detection: dumping candidates to "
             << path << endl;
    }

    std::unique_ptr<ofstream> pairDumpStream;
    std::atomic<uint64_t> pairDumpCount{0};
    if(const char* path = std::getenv("DINARA_SNP_PAIR_DUMP")) {
        pairDumpStream = std::make_unique<ofstream>(path);
    }

    // DINARA_SNP_SITE_DUMP: one line per site that survived EVERY filter,
    // "readName readPosition domAllele altAllele domCount altCount", with the
    // position in the owner read's own forward coordinates. This is what a
    // ground-truth comparison needs: the sites can be projected to reference
    // coordinates through an independent read alignment and matched against the
    // published HG002 mat-vs-pat het track. The candidate dumps above are the
    // pre-filter set and cannot measure precision or recall.
    // DINARA_SNP_OWNED_DUMP: one line per OWNED site (i.e. every site that
    // reached the filter stack), with each filter's verdict as a separate
    // column. Projected against ground truth this attributes a missed real
    // variant to the specific filter that discarded it -- which the
    // surviving-site dump cannot do, since a missed site simply is not in it.
    std::unique_ptr<ofstream> ownedDumpFile;
    ofstream* ownedDumpStream = nullptr;
    std::mutex ownedDumpMutex;
    if(const char* path = std::getenv("DINARA_SNP_OWNED_DUMP")) {
        ownedDumpFile = std::make_unique<ofstream>(path);
        (*ownedDumpFile) << "readName\tposition\tdomBase\taltBase\tdomCount"
                            "\taltCount\tsiteTotal\tbinomialPass\tccPass"
                            "\tbiallelicClean\thomopolymer\tstr\tstrandBias\n";
        ownedDumpStream = ownedDumpFile.get();
    }
    std::unique_ptr<ofstream> siteDumpFile;
    ofstream* siteDumpStream = nullptr;
    std::mutex siteDumpMutex;
    if(const char* path = std::getenv("DINARA_SNP_SITE_DUMP")) {
        siteDumpFile = std::make_unique<ofstream>(path);
        siteDumpStream = siteDumpFile.get();
    }

    std::atomic<ReadId> nextReadId{0};
    auto worker = [&]() {
        vector<uint16_t> disagree;
        vector<int32_t> coverDelta;
        vector<uint32_t> localCandidates, ownedPositions;
        vector<uint8_t> contextBuffer;
        // Per owned position, the members carrying each base. Recorded rather
        // than only counted, because these ARE the anchor arms -- an allele's
        // member list is what appendHetAnchorPair takes.
        vector<std::array<vector<pair<OrientedReadId, uint32_t>>, 4>> alleleMembers;
        vector<CigarSnpSite> localSitesOut;
        vector<uint32_t> localSitePositions;
        // Allele detail for DINARA_SNP_SITE_DUMP, index-parallel to
        // localSitePositions so the adjacency loop can dump exactly the
        // sites it keeps.
        struct SiteDetail { uint8_t domBase, altBase; uint32_t domCount, altCount; };
        vector<SiteDetail> localSiteDetail;
        // [base][strand]: strand 0 = partner aligned same-strand as the
        // owning read, 1 = reverse. A real variant should appear on both;
        // an allele seen only one way is the signature of a strand-specific
        // systematic error, which is what the Fisher test below tests for.
        vector<std::array<std::array<uint16_t, 2>, 4>> alleleCounts;
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

                // Per-record parity dump (DINARA_SNP_PAIR_DUMP): one line per
                // mismatch, "qname tname strand queryPos", for a bounded sample
                // of records. Comparing these against the same pair's line in
                // hifiasm's own PAF isolates the walk from any difference in
                // which overlaps each side has.
                if(pairDumpStream && selfIsQuery &&
                   pairDumpCount.fetch_add(1, std::memory_order_relaxed) < 200) {
                    std::ostringstream buffer;
                    const auto qn = reads->getReadName(ReadId(rec->readIdQ));
                    const auto tn = reads->getReadName(ReadId(rec->readIdT));
                    // Emit the record's whole token stream as a CIGAR string,
                    // so it can be diffed against hifiasm's own PAF line.
                    buffer.write(&*qn.begin(), std::streamsize(qn.size()));
                    buffer << '\t';
                    buffer.write(&*tn.begin(), std::streamsize(tn.size()));
                    buffer << '\t' << (rec->isSameStrand ? '+' : '-')
                           << '\t' << rec->qStart << '\t' << rec->qEnd
                           << '\t' << rec->tStart << '\t' << rec->tEnd << '\t';
                    static const char opChar[4] = {'=', 'X', 'I', 'D'};
                    for(const CigarToken tk: hifiasmImportedCigarStore.tokensOf(*rec)) {
                        buffer << tk.len() << opChar[tk.op() & 3];
                    }
                    buffer << '\n';
                    const string text = buffer.str();
                    if(!text.empty()) {
                        std::lock_guard<std::mutex> lock(dumpMutex);
                        (*pairDumpStream) << text;
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
            localCandidates.clear();
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
                    localCandidates.push_back(position);
                }
            }
            // Ownership. Every read covering a locus registers it, so the same
            // site is found once per covering read (~coverage times over).
            // Rather than emitting all of them and merging afterwards, each read
            // asks cheaply whether it OWNS the site -- lowest ReadId among the
            // reads covering that position -- and only the owner does pass 2.
            // That removes the redundancy before any expensive work, not after.
            // Cost is candidates x records interval tests, both small.
            //
            // Two reads viewing the same locus agree on the owner because they
            // minimise over the same covering set. They can disagree only where
            // their overlap sets differ.
            //
            // That disagreement is NOT only a duplicate risk: it can DROP a
            // locus. Ownership is elected over the covering reads without
            // asking whether the position is a candidate on the winner, so a
            // read that did flag the site yields to a lower ReadId that never
            // did, and then nobody runs pass 2 on it.
            //
            // Measured, that turns out to be load-bearing rather than broken,
            // and the naive fix would be a regression. On the 989-read fixture
            // 35905 distinct loci are flagged by at least one read but only
            // 3711 are ever evaluated -- ownership silently discards 90% of
            // them. Almost all of that is noise: 84.3% of the discarded loci
            // were flagged by exactly ONE read, and of the plausible real hets
            // among them there are five. Real variants are flagged
            // redundantly -- 93.8% of truth het SNVs are flagged by 10 or more
            // reads, 96.7% by at least five -- so the lowest-ReadId coverer is
            // almost always one of the flaggers and the site is kept.
            //
            // So min-ReadId election is acting as an accidental REDUNDANCY
            // filter, and it is the operative filter here: the ">= 2 disagreeing
            // partners" rule admits ten times more loci than survive this step.
            // Electing among candidate-holders instead -- the "correct" fix --
            // would admit all 32194 discarded loci into pass 2, a 9.7x increase
            // in evaluated loci that is almost entirely noise.
            //
            // The honest repair is therefore NOT to fix the election but to make
            // the redundancy requirement explicit: require a locus to be flagged
            // by >= N reads, which is what is really happening, then ownership
            // becomes pure deduplication and can be elected among flaggers
            // safely. Two things recommend it over the status quo -- the
            // threshold becomes a stated, tunable decision instead of a
            // side effect, and the result stops depending on ReadId ordering
            // (ReadIds come from input order, so which marginal loci survive
            // currently depends on the order the reads were read in).
            ownedPositions.clear();
            for(const uint32_t position: localCandidates) {
                ReadId owner = readId;
                for(const HifiasmImportedCigarStore::Record* rec: records) {
                    const bool selfIsQuery = (rec->readIdQ == readId);
                    const uint32_t from = selfIsQuery ? rec->qStart : rec->tStart;
                    const uint32_t to   = selfIsQuery ? rec->qEnd   : rec->tEnd;
                    if(position >= from && position < to) {
                        const ReadId partner =
                            selfIsQuery ? ReadId(rec->readIdT) : ReadId(rec->readIdQ);
                        if(partner < owner) owner = partner;
                    }
                }
                if(owner == readId) ownedPositions.push_back(position);
            }
            ownedSites.fetch_add(ownedPositions.size(), std::memory_order_relaxed);

            // Parity dump (DINARA_SNP_DUMP=<path>): one line per candidate,
            // "readName position nDisagree nCovering", so the counts can be
            // checked against an independent implementation over hifiasm's own
            // PAF. This is the pre-ownership candidate set, which is what is
            // comparable -- hifiasm deduplicates nothing, it corrects per read.
            if(dumpStream) {
                std::ostringstream buffer;
                for(const uint32_t position: localCandidates) {
                    int32_t cov = 0;
                    for(uint32_t q = 0; q <= position; q++) cov += coverDelta[q];
                    const auto name = reads->getReadName(readId);
                    buffer.write(&*name.begin(), std::streamsize(name.size()));
                    buffer << '\t' << position
                           << '\t' << disagree[position] << '\t' << cov << '\n';
                }
                const string text = buffer.str();
                if(!text.empty()) {
                    std::lock_guard<std::mutex> lock(dumpMutex);
                    (*dumpStream) << text;
                }
            }

            // PASS 2, for owned sites only: which base does each covering read
            // carry? The agreeing reads are one of the two allele arms, so they
            // matter as much as the disagreeing ones -- this is hifiasm's
            // addSNPtohaplotype step.
            //
            // Still no per-base loop over match runs. A run spans a contiguous
            // range of self positions, and the owned positions are sorted, so
            // the few that fall inside a run are found by binary search and
            // their partner offsets computed directly.
            if(!ownedPositions.empty()) {
                alleleCounts.assign(ownedPositions.size(), {});
                if(sitesOut) {
                    alleleMembers.assign(ownedPositions.size(), {});
                }
                // Seed with this read's own base -- it is a member too.
                for(size_t k = 0; k < ownedPositions.size(); k++) {
                    const Base self =
                        reads->getOrientedReadBase(OrientedReadId(readId, 0), ownedPositions[k]);
                    if(self.value < 4) {
                        alleleCounts[k][self.value][0]++;
                        if(sitesOut) {
                            alleleMembers[k][self.value].emplace_back(
                                OrientedReadId(readId, 0), ownedPositions[k]);
                        }
                    }
                }

                for(const HifiasmImportedCigarStore::Record* rec: records) {
                    const bool selfIsQuery = (rec->readIdQ == readId);
                    const uint32_t targetLength =
                        uint32_t(reads->getRead(ReadId(rec->readIdT)).baseCount);
                    // The partner, in the orientation the alignment sees it.
                    // The query side is always forward; the target side is
                    // reverse-complemented exactly when the overlap is.
                    const OrientedReadId partner = selfIsQuery ?
                        OrientedReadId(ReadId(rec->readIdT), rec->isSameStrand ? 0 : 1) :
                        OrientedReadId(ReadId(rec->readIdQ), 0);
                    uint64_t qPos = rec->qStart;
                    uint64_t tPos = rec->isSameStrand ?
                        rec->tStart : (targetLength - rec->tEnd);

                    for(const CigarToken token: hifiasmImportedCigarStore.tokensOf(*rec)) {
                        const uint8_t op = token.op();
                        const uint16_t length = token.len();
                        if(op == CigarOpMatch || op == CigarOpMismatch) {
                            // Self positions covered by this run, as a range.
                            // Reverse-strand targets run backwards.
                            const bool descending = (!selfIsQuery && !rec->isSameStrand);
                            const int64_t first = selfIsQuery ? int64_t(qPos) :
                                (rec->isSameStrand ? int64_t(tPos) :
                                 int64_t(targetLength) - 1 - int64_t(tPos));
                            const int64_t last = descending ?
                                (first - int64_t(length) + 1) : (first + int64_t(length) - 1);
                            const int64_t lo = std::min(first, last);
                            const int64_t hi = std::max(first, last);

                            auto it = std::lower_bound(ownedPositions.begin(),
                                ownedPositions.end(), uint32_t(std::max<int64_t>(0, lo)));
                            for(; it != ownedPositions.end() && int64_t(*it) <= hi; ++it) {
                                const int64_t offset = descending ?
                                    (first - int64_t(*it)) : (int64_t(*it) - first);
                                if(offset < 0 || offset >= int64_t(length)) continue;
                                const uint32_t partnerPosition = uint32_t(
                                    (selfIsQuery ? tPos : qPos) + uint64_t(offset));
                                Base base =
                                    reads->getOrientedReadBase(partner, partnerPosition);
                                // Express the partner's base in THIS read's
                                // forward frame. When this read is the target of
                                // a reverse-strand overlap the alignment sees it
                                // reverse-complemented, so an alignment match
                                // means partner == complement(our forward base).
                                // Without this the allele identity is the
                                // complement of what it should be.
                                if(!selfIsQuery && !rec->isSameStrand && base.value < 4) {
                                    base.value = uint8_t(3 - base.value);
                                }
                                if(base.value < 4) {
                                    const size_t slot = size_t(it - ownedPositions.begin());
                                    alleleCounts[slot][base.value]
                                        [rec->isSameStrand ? 0 : 1]++;
                                    if(sitesOut) {
                                        alleleMembers[slot][base.value].emplace_back(
                                            partner, partnerPosition);
                                    }
                                }
                                // Self-check against the reads themselves, which
                                // needs no second hifiasm run and so cannot be
                                // confounded by seeding differences. At a MATCH
                                // column the two bases must be identical; at a
                                // MISMATCH column they must differ. The first
                                // catches bad orientation handling, the second
                                // catches a mismatch we report that hifiasm's
                                // CIGAR and the sequence do not support -- which
                                // is exactly what "full parity for SNPs" means.
                                const Base self = reads->getOrientedReadBase(
                                    OrientedReadId(readId, 0), *it);
                                if(op == CigarOpMatch) {
                                    matchColumnsChecked.fetch_add(1, std::memory_order_relaxed);
                                    if(self.value == base.value) {
                                        matchColumnsAgree.fetch_add(1, std::memory_order_relaxed);
                                    }
                                } else {
                                    mismatchColumnsChecked.fetch_add(1, std::memory_order_relaxed);
                                    if(self.value != base.value) {
                                        mismatchColumnsDiffer.fetch_add(1, std::memory_order_relaxed);
                                    }
                                }
                            }
                        }
                        if(opConsumesQuery(op))  qPos += length;
                        if(opConsumesTarget(op)) tPos += length;
                    }
                }

                // How many distinct bases does each site actually carry, and
                // does the site survive the significance test?
                //
                // The raw ">= 2 reads" count is a fixed floor and a poor one:
                // it accepts two independent errors that happen to agree. The
                // binomial test asks the right question instead -- given this
                // many reads on the dominant allele and an assumed per-read
                // error rate, how surprising is it that k reads carry the
                // minority one? That is the same test the abPOA detector
                // applies (now shared, see hetSignificance.hpp), so both
                // detectors agree on what counts as a site.
                for(size_t k = 0; k < alleleCounts.size(); k++) {
                    const auto& strandCounts = alleleCounts[k];
                    std::array<uint16_t, 4> counts{};
                    for(int b = 0; b < 4; b++) {
                        counts[b] = uint16_t(strandCounts[b][0] + strandCounts[b][1]);
                    }
                    uint64_t distinct = 0, total = 0, dominant = 0;
                    for(const uint16_t c: counts) {
                        total += c;
                        if(c >= 2) distinct++;   // a lone read is not an allele
                        if(c > dominant) dominant = c;
                    }

                    // Sequence-context gates, evaluated on THIS read around the
                    // site. Statistics cannot save us here: ONT homopolymer and
                    // short-tandem-repeat errors are systematic, so reads
                    // disagree there consistently and the site looks
                    // convincingly biallelic to any significance test. Only
                    // context can reject it. Same primitive the window path
                    // uses -- unit length 1 is a homopolymer, 2..6 an STR.
                    const uint32_t sitePosition = ownedPositions[k];
                    bool inHomopolymer = false, inStr = false;
                    {
                        const uint32_t lo = (sitePosition > 16) ? (sitePosition - 16) : 0;
                        const uint32_t hi = std::min<uint32_t>(readLength, sitePosition + 17);
                        contextBuffer.clear();
                        for(uint32_t q = lo; q < hi; q++) {
                            contextBuffer.push_back(
                                reads->getOrientedReadBase(OrientedReadId(readId, 0), q).value);
                        }
                        KmVarKey key{};
                        key.type = KmVarType::Snp;
                        key.pos = uint32_t(sitePosition - lo);
                        inHomopolymer = kmIsRepeatUnitRange(
                            contextBuffer.data(), uint32_t(contextBuffer.size()), key, 0, 1, 1);
                        inStr = kmIsRepeatUnitRange(
                            contextBuffer.data(), uint32_t(contextBuffer.size()), key, 0, 2, 6);
                    }
                    // The homopolymer gate is measured to be NET HARMFUL and is
                    // off by default. Against the HG002 truth track it was the
                    // single largest cause of missed real variants -- 49 of the
                    // 58 truth SNVs that reached the filters and were rejected
                    // died here -- while removing almost no false positives,
                    // because the VAF floor above already accounts for those.
                    // Disabling it takes recall 57.8% -> 70.0% with precision
                    // unchanged (97.5% -> 98.0%), on the naive all-truth
                    // denominator in use at the time. Note the circularity this
                    // avoids: scoring a homopolymer gate against truth calls
                    // that are themselves unreliable in homopolymers measures
                    // agreement between two guesses, which is why the truth
                    // partition at the top of this file excludes them.
                    //
                    // The rationale for having it was sound (ONT homopolymer
                    // error is systematic, so no significance test can reject
                    // it) but incomplete: real heterozygous SNVs are common
                    // immediately beside homopolymers, and a +/-16 bp context
                    // window rejects those too. The STR gate (unit 2..6) stays
                    // on: dropping it as well adds recall (-> 73.6%) but costs
                    // precision (-> 94.4%), a bad trade when each surviving site
                    // becomes an anchor.
                    if(inHomopolymer && !filterHomopolymer) inHomopolymer = false;
                    if(inStr && !filterStr) inStr = false;
                    if(inHomopolymer) droppedHomopolymer.fetch_add(1, std::memory_order_relaxed);
                    if(inStr)         droppedRepeat.fetch_add(1, std::memory_order_relaxed);
                    allelesTotal.fetch_add(total, std::memory_order_relaxed);
                    alleleCountHistogram[std::min<size_t>(8, size_t(distinct))]
                        .fetch_add(1, std::memory_order_relaxed);

                    // The dominant allele is the implicit reference and is
                    // never itself tested; every other allele must clear the
                    // tail probability.
                    uint64_t passing = 0;
                    bool dominantSeen = false;
                    for(const uint16_t c: counts) {
                        if(c == 0) continue;
                        if(!dominantSeen && c == dominant) { dominantSeen = true; ++passing; continue; }
                        if(binomialTailPValue(dominant, c, hetErrorRate) <= hetSignificance) {
                            ++passing;
                        }
                    }
                    if(passing >= 2) significantSites.fetch_add(1, std::memory_order_relaxed);
                    significanceHistogram[std::min<size_t>(4, size_t(passing))]
                        .fetch_add(1, std::memory_order_relaxed);

                    // Reduce the site to biallelic, or reject it -- hifiasm's
                    // rule (Correct.cpp, the gates before InsertSNPVector). It
                    // does not split a multi-allelic site into several binary
                    // ones; it picks the single strongest alternate and demands
                    // the site be effectively biallelic, dropping it otherwise.
                    // The binomial test above says an allele is not noise; these
                    // say the SITE is clean enough to phase on.
                    uint64_t bestAlt = 0; int bestAltBase = -1; uint64_t ties = 0;
                    {
                        bool dominantTaken = false;
                        for(int b = 0; b < 4; b++) {
                            if(counts[b] == 0) continue;
                            if(!dominantTaken && counts[b] == dominant) { dominantTaken = true; continue; }
                            if(counts[b] > bestAlt) { bestAlt = counts[b]; bestAltBase = b; }
                        }
                        // Two alternates tied at the top: which one is the other
                        // haplotype is ambiguous, so hifiasm drops the site.
                        if(bestAltBase >= 0) {
                            bool takenAgain = false;
                            for(int b = 0; b < 4; b++) {
                                if(counts[b] == 0) continue;
                                if(!takenAgain && counts[b] == dominant) { takenAgain = true; continue; }
                                if(b != bestAltBase && counts[b] == bestAlt) ties++;
                            }
                        }
                    }
                    bool biallelicClean = (bestAltBase >= 0) && (bestAlt >= 2) && (ties == 0);
                    if(biallelicClean) {
                        // Reference plus the chosen alternate must account for
                        // essentially the whole pileup.
                        if(double(dominant + bestAlt) / double(total > 0 ? total : 1)
                           < siteMinPurity) {
                            biallelicClean = false;
                            droppedImpure.fetch_add(1, std::memory_order_relaxed);
                        } else {
                            // And the chosen alternate must dominate the
                            // disagreement, not merely lead a scattered field.
                            const uint64_t disagreeing = (total > dominant) ? (total - dominant) : 0;
                            if(disagreeing > 0 &&
                               double(bestAlt) / double(disagreeing) < siteMinAltDominance) {
                                biallelicClean = false;
                                droppedNotDominant.fetch_add(1, std::memory_order_relaxed);
                            }
                        }
                    } else if(ties > 0) {
                        droppedTiedAlleles.fetch_add(1, std::memory_order_relaxed);
                    }

                    // Strand bias, on the strongest minority allele. A real
                    // variant is seen from both directions; one that is not is
                    // the signature of a strand-specific systematic error.
                    // Same form as the longcallD-derived test already in the
                    // tree: Fisher exact on (fwdAlt, revAlt) against a balanced
                    // expectation.
                    bool strandBiased = false;
                    if(passing >= 2) {
                        uint64_t altBest = 0; int altBase = -1;
                        bool dominantTaken = false;
                        for(int b = 0; b < 4; b++) {
                            if(counts[b] == 0) continue;
                            if(!dominantTaken && counts[b] == dominant) { dominantTaken = true; continue; }
                            if(counts[b] > altBest) { altBest = counts[b]; altBase = b; }
                        }
                        if(altBase >= 0) {
                            const int fa = int(strandCounts[altBase][0]);
                            const int ra = int(strandCounts[altBase][1]);
                            const int expected = (fa + ra) / 2;
                            if(expected > 0) {
                                const double p = kmFisherExactTwoTail(fa, ra, expected, expected);
                                if(p < strandBiasPValue) {
                                    strandBiased = true;
                                    droppedStrandBias.fetch_add(1, std::memory_order_relaxed);
                                }
                            }
                        }
                    }

                    // Absolute depth floor on the MAJOR allele: hifiasm's
                    // cc, and nothing else.
                    //
                    // The binomial test above and this floor interrogate
                    // DIFFERENT alleles, which is why both are needed and why
                    // neither subsumes the other. The binomial asks whether the
                    // MINOR allele is too large to be a misread of the dominant
                    // one -- a purely RELATIVE question, with no notion of
                    // absolute depth. That is exactly where it fails: at a
                    // dominant count of 1, P(X >= 1) at the assumed error rate
                    // IS that rate, so a 1-vs-1 site clears p <= 0.05 wherever
                    // coverage is thin. cc asks the complementary question --
                    // does the MAJOR allele carry the support a real haplotype
                    // would, measured against this dataset's own het coverage
                    // peak -- and so supplies precisely the depth the binomial
                    // test cannot see.
                    //
                    // This replaces a pair of VAF thresholds previously applied
                    // here, whose shape came from hifiasm's
                    // filter_one_snp_advance_nearby. Those are gone deliberately:
                    // that function is part of hifiasm's LEGACY clustering path
                    // (reachable only from cluster_advance / cluster_ul_advance,
                    // never from the live rphase_hc EC path), it filters whole
                    // candidate haplotype VECTORS inside a path DP rather than
                    // one site, and its constants were tuned against a
                    // fraction-window prefilter this detector does not have. cc
                    // is from the live path, is derived from the run's measured
                    // coverage rather than fixed, and gates the quantity the
                    // binomial test leaves ungated -- so layering the VAF
                    // clauses on top of it only removed sites twice for the same
                    // reason, with the weaker justification.
                    // A minor-allele FRACTION floor sits beside cc. Unlike the
                    // two-tier VAF rule this replaced, it is justified by
                    // measurement rather than by hifiasm provenance: projected
                    // against the published HG002 mat-vs-pat het track on the
                    // 989-read fixture, the false positives pile up at a VAF
                    // near 0.1 (48 of 58) while real het variants sit at 0.4-0.5
                    // (211 of 236), so the two populations barely overlap. A
                    // floor of 0.25 removes essentially all of the former at
                    // almost no cost to the latter: precision 75.4% -> 97.5%
                    // while recall was flat (58.7% -> 57.8% on the naive
                    // all-truth denominator in use at the time; see the
                    // MEASURED ACCURACY note at the top of this file for the
                    // corrected denominators and the flat threshold surface).
                    //
                    // These low-VAF sites clear both other tests legitimately --
                    // 40-vs-4 passes cc easily and is significant under the
                    // binomial at an 0.025 error rate -- so neither the depth
                    // floor nor the significance test can stand in for this one.
                    // They are recurrent systematic error or paralog signal, and
                    // only the allele BALANCE separates them from real hets.
                    bool enoughCoverage = false, ccPass = false, vafPass = false;
                    if(bestAltBase >= 0) {
                        ccPass = (dominant >= minAlleleCoverage);
                        const uint64_t pairTotal = dominant + bestAlt;
                        const double vaf = (pairTotal > 0) ?
                            (double(bestAlt) / double(pairTotal)) : 0.0;
                        vafPass = (vaf >= minAlleleFraction);
                        enoughCoverage = ccPass && vafPass;
                    }
                    if(!enoughCoverage) {
                        droppedLowCoverage.fetch_add(1, std::memory_order_relaxed);
                    }

                    // How much do the binomial test and cc actually overlap?
                    // They interrogate different alleles -- the binomial asks
                    // whether the MINOR allele is too large to be error, cc
                    // whether the MAJOR allele is large enough to be a
                    // haplotype -- so they should reject largely different
                    // sites. Counted here over every owned site, independent of
                    // the other gates, purely to see that.
                    // Restricted to sites that HAVE an alternate allele: a
                    // site with none has no minor allele for either test to
                    // interrogate, and counting it would score as a cc
                    // rejection something neither test is really deciding.
                    if(bestAltBase >= 0) {
                        const bool binomialPass = (passing >= 2);
                        const size_t cell = (binomialPass ? 2 : 0) + (ccPass ? 1 : 0);
                        binomialVsCc[cell].fetch_add(1, std::memory_order_relaxed);
                    } else {
                        noAlternateSites.fetch_add(1, std::memory_order_relaxed);
                    }

                    if(ownedDumpStream) {
                        uint64_t siteTotal = 0;
                        for(int b = 0; b < 4; b++) siteTotal += counts[b];
                        int domB = -1;
                        for(int b = 0; b < 4; b++) if(counts[b] == dominant) { domB = b; break; }
                        const auto nm = reads->getReadName(readId);
                        std::ostringstream line;
                        line.write(&*nm.begin(), std::streamsize(nm.size()));
                        line << '\t' << sitePosition
                             << '\t' << (domB >= 0 ? Base::fromInteger(uint8_t(domB)).character() : '.')
                             << '\t' << (bestAltBase >= 0 ? Base::fromInteger(uint8_t(bestAltBase)).character() : '.')
                             << '\t' << dominant << '\t' << bestAlt << '\t' << siteTotal
                             << '\t' << (passing >= 2 ? 1 : 0)
                             << '\t' << (ccPass ? 1 : 0)
                             << '\t' << (biallelicClean ? 1 : 0)
                             << '\t' << (inHomopolymer ? 1 : 0)
                             << '\t' << (inStr ? 1 : 0)
                             << '\t' << (strandBiased ? 1 : 0)
                             << '\n';
                        const string text = line.str();
                        std::lock_guard<std::mutex> lock(ownedDumpMutex);
                        (*ownedDumpStream) << text;
                    }

                    if(passing >= 2 && biallelicClean && enoughCoverage &&
                       !inHomopolymer && !inStr && !strandBiased) {
                        sitesAfterFilters.fetch_add(1, std::memory_order_relaxed);
                        if(sitesOut) {
                            // Two arms: the dominant allele and the chosen
                            // alternate. Reads in neither get no new anchor and
                            // simply take the direct flank-to-flank edge once
                            // journeys are rebuilt -- the same contract the
                            // abPOA detector relies on.
                            int dominantBase = -1;
                            for(int b = 0; b < 4; b++) {
                                if(counts[b] == dominant) { dominantBase = b; break; }
                            }
                            if(dominantBase >= 0 && bestAltBase >= 0) {
                                CigarSnpSite site;
                                site.alleles.push_back(alleleMembers[k][dominantBase]);
                                site.alleles.push_back(alleleMembers[k][bestAltBase]);
                                localSitesOut.push_back(std::move(site));
                                localSitePositions.push_back(sitePosition);
                                localSiteDetail.push_back({uint8_t(dominantBase),
                                    uint8_t(bestAltBase), uint32_t(dominant),
                                    uint32_t(bestAlt)});
                            }
                        }
                    }
                }
            }

            if(sitesOut && !localSitesOut.empty()) {
                // Adjacent-site rejection, from the LIVE ONT path
                // (generate_haplotypes_naive_HiFi's "filter snps" loop, which
                // drops any site whose neighbour sits at site +/- 1). Two het
                // sites at consecutive bases are almost always one misplaced
                // indel spelled as two substitutions -- exactly the
                // substitution-vs-indel spelling ambiguity the CIGAR
                // comparison turned up -- rather than two independent
                // variants. Positions here are all on this one read and were
                // produced in ascending order, so neighbours are adjacent
                // entries.
                for(size_t k = 0; k < localSitePositions.size(); k++) {
                    const bool prevAdjacent = (k > 0) &&
                        (localSitePositions[k] == localSitePositions[k - 1] + 1);
                    const bool nextAdjacent = (k + 1 < localSitePositions.size()) &&
                        (localSitePositions[k] + 1 == localSitePositions[k + 1]);
                    if(prevAdjacent || nextAdjacent) {
                        droppedAdjacent.fetch_add(1, std::memory_order_relaxed);
                        continue;
                    }
                    if(siteDumpStream) {
                        const auto name = reads->getReadName(readId);
                        const SiteDetail& d = localSiteDetail[k];
                        std::ostringstream line;
                        line.write(&*name.begin(), std::streamsize(name.size()));
                        line << '\t' << localSitePositions[k]
                             << '\t' << Base::fromInteger(d.domBase)
                             << '\t' << Base::fromInteger(d.altBase)
                             << '\t' << d.domCount << '\t' << d.altCount << '\n';
                        const string text = line.str();
                        std::lock_guard<std::mutex> lock(siteDumpMutex);
                        (*siteDumpStream) << text;
                    }
                    std::lock_guard<std::mutex> lock(sitesMutex);
                    sitesOut->push_back(std::move(localSitesOut[k]));
                }
                localSitesOut.clear();
                localSitePositions.clear();
                localSiteDetail.clear();
            } else if(sitesOut) {
                localSitePositions.clear();
                localSiteDetail.clear();
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
    cout << "  distinct sites owned (deduplicated across covering reads): "
         << ownedSites.load() << endl;
    {
        const uint64_t checked = matchColumnsChecked.load();
        const uint64_t agree = matchColumnsAgree.load();
        const uint64_t mchecked = mismatchColumnsChecked.load();
        const uint64_t mdiffer = mismatchColumnsDiffer.load();
        cout << "  mismatch-column self-check: " << mdiffer << " / " << mchecked
             << " really differ ("
             << (mchecked ? 100.0*double(mdiffer)/double(mchecked) : 0.0)
             << "%) -- must be 100%, else we report SNPs the sequence "
                "does not support" << endl;
        cout << "  match-column self-check: " << agree << " / " << checked
             << " agree (" << (checked ? 100.0*double(agree)/double(checked) : 0.0)
             << "%) -- must be 100%, else the orientation handling is wrong" << endl;
    }
    cout << "  sites passing the binomial significance test (>= 2 alleles "
            "significant at p <= " << hetSignificance << ", error rate "
         << hetErrorRate << "): " << significantSites.load() << " of "
         << ownedSites.load() << endl;
    cout << "  significant alleles per owned site:" << endl;
    for(size_t i = 0; i < significanceHistogram.size(); i++) {
        const uint64_t count = significanceHistogram[i].load();
        if(count) cout << "    " << i << ": " << count << endl;
    }
    cout << "  context/strand filters: " << droppedHomopolymer.load()
         << " in a homopolymer, " << droppedRepeat.load() << " in an STR, "
         << droppedStrandBias.load() << " strand-biased (Fisher p < "
         << strandBiasPValue << ")" << endl;
    cout << "  biallelic reduction (hifiasm's rule): "
         << droppedTiedAlleles.load() << " with tied alternates, "
         << droppedImpure.load() << " below " << siteMinPurity << " purity, "
         << droppedNotDominant.load() << " alternate below "
         << siteMinAltDominance << " of the disagreement" << endl;
    cout << "  below the major-allele depth floor (cc = " << minAlleleCoverage
         << ") or the minor-allele fraction floor (" << minAlleleFraction
         << "): " << droppedLowCoverage.load() << endl;
    cout << "  binomial vs cc, over all owned sites (they gate different "
            "alleles, so they should reject different sites):"
         << "\n    rejected by both          : " << binomialVsCc[0].load()
         << "\n    rejected by binomial only : " << binomialVsCc[1].load()
         << "\n    rejected by cc only       : " << binomialVsCc[2].load()
         << "\n    passed by both            : " << binomialVsCc[3].load()
         << "\n    (excludes " << noAlternateSites.load()
         << " owned sites with no alternate allele at all)" << endl;
    cout << "  adjacent to another site (+/- 1 bp, hifiasm's live ONT rule): "
         << droppedAdjacent.load() << endl;
    cout << "  sites surviving ALL filters: " << sitesAfterFilters.load()
         << " of " << ownedSites.load() << endl;
    cout << "  alleles per owned site (a base needs >= 2 reads to count):" << endl;
    for(size_t i = 0; i < alleleCountHistogram.size(); i++) {
        const uint64_t count = alleleCountHistogram[i].load();
        if(count) cout << "    " << i << " alleles: " << count << endl;
    }
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
