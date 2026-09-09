#include "Shasta2Journeys.hpp"
#include "Shasta2Anchors.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

#include <algorithm>
#include <iostream>
#include <thread>
using std::cout;
using std::endl;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
namespace dinara {
    template class MultithreadedObject<Shasta2Journeys>;
}



// Initial creation.
// This sets the positionInJourney for every AnchorMarkerInfo
// stored in the Anchors, and for this reason the Anchors
// are not passed in as const.
Shasta2Journeys::Shasta2Journeys(
    uint64_t orientedReadCount,
    shared_ptr<Shasta2Anchors> anchorsPointer,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    anchorsPointer(anchorsPointer)

{
    performanceLog << timestamp << "Journeys creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t anchorCount = anchorsPointer->size();
    const uint64_t anchorBatchCount = 1000;
    const uint64_t orientedReadBatchCount = 1000;

    // Pass1: make space for the journeysWithOrdinals.
    journeysWithOrdinals.createNew(largeDataName("tmp-Shasta2JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Shasta2Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Shasta2Journeys::threadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithOrdinals to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Shasta2Journeys::threadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "Journeys creation ends." << endl;
}



void Shasta2Journeys::threadFunction1(uint64_t /* threadId */)
{
    threadFunction12(1);
}



void Shasta2Journeys::threadFunction2(uint64_t /* threadId */)
{
    threadFunction12(2);
}



void Shasta2Journeys::threadFunction12(uint64_t pass)
{
    const Shasta2Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(Shasta2AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            Shasta2Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.ordinal});
                }
            }
        }
    }
}



void Shasta2Journeys::threadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithOrdinals[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Shasta2Journeys::threadFunction4(uint64_t /* threadId */)
{
    Shasta2Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Copy the journeysWithOrdinals to the journeys.
            const auto v = journeysWithOrdinals[orientedReadValue];
            const auto journey = journeys[orientedReadValue];
            DINARA_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            // Anchor members are stored sorted ascending by OrientedReadId (see
            // Shasta2Anchors), so binary search instead of scanning the whole anchor.
            for(uint64_t position=0; position<journey.size(); position++) {
                const Shasta2AnchorId anchorId = journey[position];
                span<Shasta2AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
                const auto it = std::lower_bound(markerInfos.begin(), markerInfos.end(), orientedReadId,
                    [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                        return info.orientedReadId < oid;
                    });
                DINARA_ASSERT(it != markerInfos.end() and it->orientedReadId == orientedReadId);
                it->positionInJourney = uint32_t(position);
            }
        }
    }
}











// Rebuild the journeys and positionInJourney from scratch using the current
// anchor set, keyed by `position` instead of `ordinal` -- see the .hpp
// comment for why. Strand 1 is never sorted independently: it is derived by
// reversing strand 0's sorted list and flipping each anchor id, exactly like
// the RC mirror emitted below, so the journey(R,1) ==
// reverse(journey(R,0)) invariant holds regardless of how same-position ties
// on strand 0 happen to be broken (ascending anchorId here, arbitrary but
// deterministic). Runs serially: the anchor count here is small (tens of
// thousands), so a second full pass is cheap and avoids any thread-safety
// concerns in a rebuild that runs once, off the hot path.
void Shasta2Journeys::rebuildAfterNewAnchors(Shasta2AnchorId newAnchorsBegin, uint64_t threadCount)
{
    performanceLog << timestamp << "Journeys rebuild (new anchors) begins." << endl;
    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeys.isOpen());
    static_cast<void>(threadCount);

    Shasta2Anchors& anchors = *anchorsPointer;
    const uint64_t anchorCount = anchors.size();
    const uint64_t orientedReadCount = journeys.size();
    const uint64_t readCount = orientedReadCount / 2;
    const Shasta2AnchorId hetFirst = anchors.hetAnchorFirstId;

    // Pass A: collect, per read (strand 0 only), (position, anchorId) pairs by
    // scanning primary anchors and the new anchor range only -- see the .hpp
    // comment for why anchors in between (e.g. export-only SNPmer anchors)
    // are excluded.
    vector<vector<pair<uint32_t, Shasta2AnchorId>>> strand0(readCount);
    for(Shasta2AnchorId anchorId = 0; anchorId < anchorCount; anchorId++) {
        const bool isPrimary = (hetFirst == invalid<Shasta2AnchorId>) || (anchorId < hetFirst);
        const bool isNew = (anchorId >= newAnchorsBegin);
        if(!isPrimary && !isNew) continue;
        for(const Shasta2AnchorMarkerInfo& info : anchors.anchorMarkerInfos[anchorId]) {
            if((info.orientedReadId.getValue() & 1U) != 0U) continue;   // strand 1 done via mirror.
            strand0[info.orientedReadId.getReadId()].push_back({info.position, anchorId});
        }
    }

    // Pass B: sort each read's strand-0 list by position, then collapse any
    // EXACT position tie down to a single survivor. Two different anchors
    // legitimately CAN land on the same raw position on the same read: the
    // per-edge margin proof (see Shasta2AnchorGraphHetOnGraph.cpp's file
    // header) only guarantees separation from an edge's OWN flanking
    // anchors, not from anchors created by some OTHER, independent edge that
    // also happens to touch this read at the same base (common once real
    // data has many overlapping anchor-graph edges) -- left unresolved, two
    // such anchors land journey-adjacent with a zero offset and
    // Shasta2AnchorPair::assertNoNegativeOffsets aborts the whole run. The
    // loser is dropped from just THIS read's journey -- it keeps its other
    // members and its own coverage count untouched, exactly like the
    // export-side drop map.
    //
    // WHICH one survives is the only freedom available (shasta2's
    // LocalAssembly7::positionOffsetAB() asserts positionB > positionA
    // strictly, so keeping both is not an option), and the choice is not
    // symmetric. Preferring primary over het -- which is what this did, and
    // what main.cpp's export-side resolution still does -- makes het anchors
    // lose EVERY collision they are in, because a het anchor is by definition
    // never the primary. Measured on the 989-read fixture: 1154 of 1154 dropped
    // occurrences were het, none primary, costing 11.1% of all het-anchor
    // occurrences and diluting every site's arms by about that much.
    //
    // Preferring het inverts a systematic loss into a diffuse one. The
    // occurrence still has to be dropped, but it now lands on a primary anchor
    // carrying roughly twice the coverage, where one lost occurrence is
    // proportionally far cheaper -- and het anchors exist precisely to be
    // traversed, so dropping them defeats the reason they were created. Note
    // this changes journeys only: the primary anchor's own member count and
    // coverage are untouched, so the graph is still verified consistent with
    // the journeys either way (DINARA_VERIFY_ANCHOR_GRAPH).
    const bool preferHetOnTie = journeyTiePreferHet;

    // Occurrences still standing per anchor, so tie resolution can never strip
    // an anchor of its LAST one. Without this an anchor can be evicted from
    // every journey it appears in -- measured: exactly 4 low-coverage primary
    // anchors on the 989-read fixture, whose every occurrence collides with a
    // het anchor and therefore loses every tie under any class preference.
    //
    // Such an anchor is not merely diluted, it disappears: no journey contains
    // it, so it gets no edges and the anchor graph has a vertex nothing
    // traverses. Nothing upstream creates an empty anchor -- primary anchors
    // clear minAnchorCoverage and het arms clear their own >= 2 floor -- so an
    // anchor with no occurrences is manufactured here and nowhere else.
    //
    // The guard has to live at THIS layer rather than at the export. The export
    // resolves the same ties again on the anchor objects, and guarding only
    // there produced an anchor that had members in the exported anchor set but
    // no occurrence in any journey, so the two exported artifacts disagreed.
    vector<uint64_t> remainingOccurrences(anchorCount, 0);
    for(uint64_t readIdValue = 0; readIdValue < readCount; readIdValue++) {
        for(const auto& [position, anchorId] : strand0[readIdValue]) {
            static_cast<void>(position);
            remainingOccurrences[anchorId]++;
        }
    }
    // Only anchors that HAD an occurrence can be evicted from all of them.
    // Comparing against zero instead would count every anchor with no strand-0
    // occurrence at all (RC anchors, which are mirrored rather than listed).
    const vector<uint64_t> initialOccurrences = remainingOccurrences;
    // Members left on each CANONICAL anchor once the drops recorded so far are
    // applied. A tie drop removes one member from the canonical anchor (and its
    // RC twin), so this is what the export will actually write.
    vector<uint64_t> remainingMembers(anchorCount, 0);
    for(uint64_t a = 0; a < anchorCount; a++) remainingMembers[a] = anchors[a].size();
    auto isHetAnchor = [&](Shasta2AnchorId id) -> bool {
        return hetFirst != invalid<Shasta2AnchorId> && id >= hetFirst;
    };
    auto keepAOverB = [&](Shasta2AnchorId a, Shasta2AnchorId b) -> bool {
        // An anchor that cannot spare a member wins outright, ahead of any
        // class preference: losing here would take it below the coverage floor
        // it was selected under. Everything else is free to lose one.
        const uint64_t canonA = a & ~Shasta2AnchorId(1);
        const uint64_t canonB = b & ~Shasta2AnchorId(1);
        const bool aLast = (remainingOccurrences[a] <= 1) ||
                           (remainingMembers[canonA] <= minAnchorCoverage);
        const bool bLast = (remainingOccurrences[b] <= 1) ||
                           (remainingMembers[canonB] <= minAnchorCoverage);
        if(aLast != bLast) return aLast;
        const bool aHet = isHetAnchor(a), bHet = isHetAnchor(b);
        if(aHet != bHet) return preferHetOnTie ? aHet : !aHet;
        const uint64_t ca = anchors[a].size(), cb = anchors[b].size();
        if(ca != cb) return ca > cb;                         // higher coverage wins
        return a < b;                                        // lower id wins
    };
    uint64_t tieGroups = 0, unitsDropped = 0;
    journeyTieDrops.clear();
    // A het anchor dropped here never reaches the anchor graph for this read,
    // which defeats the reason it was created -- so count the losers by kind
    // and report the rate, because the tie-break makes the loss SYSTEMATIC
    // rather than incidental: "primary over het/hom" means a het anchor loses
    // every collision it is in. Measured on the 989-read fixture: 1154 of 1154
    // dropped occurrences were het anchors and none were primary, costing 11.1%
    // of all het-anchor occurrences (1154 lost against 9228 surviving), which
    // dilutes each site's arms by roughly that fraction.
    //
    // It is forced, not a choice that can simply be reversed here: shasta2's
    // LocalAssembly7::positionOffsetAB() asserts positionB > positionA
    // STRICTLY, so two anchors at the same base in one read cannot be
    // journey-adjacent at all. Keeping both is not available; only which one
    // survives is. Inverting the priority would trade an 11% dilution of het
    // arms for a much smaller proportional loss on primary anchors (they carry
    // ~2x the coverage), but primary anchors are the graph backbone, so that
    // swap needs its own verification run before it could be trusted.
    uint64_t hetDropped = 0, primaryDropped = 0, hetOccurrences = 0;
    // Why a het occurrence lost: the only way it can, given het is preferred,
    // is that the anchor beating it could not spare a member without falling
    // below its coverage floor.
    uint64_t hetLostToFloor = 0, hetLostOther = 0;
    // Original member count of the anchor that beat a het occurrence, to tell a
    // genuinely marginal anchor (selected at coverage 2) from one whittled down
    // to the floor by earlier drops -- the second would make the outcome depend
    // on the order ties are visited.
    vector<uint64_t> floorWinnerOriginalCoverage;
    // Size of the tie group, and how many of its members were het.
    vector<uint64_t> floorGroupSize, floorGroupHets;
    vector<uint64_t> overlapShared, overlapPrimary;

    filteredJourneys.assign(orientedReadCount, {});
    for(uint64_t readIdValue = 0; readIdValue < readCount; readIdValue++) {
        auto& v = strand0[readIdValue];
        std::sort(v.begin(), v.end());

        if(!v.empty()) {
            vector<pair<uint32_t, Shasta2AnchorId>> deduped;
            deduped.reserve(v.size());
            size_t i = 0;
            while(i < v.size()) {
                size_t j = i + 1;
                while(j < v.size() && v[j].first == v[i].first) j++;
                if(j - i == 1) {
                    deduped.push_back(v[i]);
                } else {
                    ++tieGroups;
                    Shasta2AnchorId keeper = v[i].second;
                    for(size_t t = i + 1; t < j; t++) {
                        if(keepAOverB(v[t].second, keeper)) keeper = v[t].second;
                    }
                    unitsDropped += (j - i - 1);
                    for(size_t t = i; t < j; t++) {
                        if(v[t].second == keeper) continue;
                        if(isHetAnchor(v[t].second)) {
                            hetDropped++;
                            const uint64_t kc = keeper & ~Shasta2AnchorId(1);
                            if(!isHetAnchor(keeper) &&
                               (remainingMembers[kc] <= minAnchorCoverage ||
                                remainingOccurrences[keeper] <= 1)) {
                                hetLostToFloor++;
                                floorWinnerOriginalCoverage.push_back(anchors[kc].size());
                                floorGroupSize.push_back(j - i);
                                uint64_t nHet = 0;
                                for(size_t u = i; u < j; u++)
                                    if(isHetAnchor(v[u].second)) nHet++;
                                floorGroupHets.push_back(nHet);
                                // How much of the primary anchor's membership
                                // does the het anchor that displaced it share?
                                // If it is nearly all of it, the two describe
                                // the same locus and the het arm is redundant
                                // with an anchor that already exists.
                                {
                                    const uint64_t hc = v[t].second & ~Shasta2AnchorId(1);
                                    uint64_t shared = 0;
                                    for(const Shasta2AnchorMarkerInfo& mp : anchors[kc])
                                        for(const Shasta2AnchorMarkerInfo& mh : anchors[hc])
                                            if(mp.orientedReadId == mh.orientedReadId) { shared++; break; }
                                    overlapShared.push_back(shared);
                                    overlapPrimary.push_back(anchors[kc].size());
                                }
                            } else {
                                hetLostOther++;
                            }
                        }
                        else primaryDropped++;
                        if(remainingOccurrences[v[t].second] > 0)
                            remainingOccurrences[v[t].second]--;
                        // The unit the export must omit. Canonical (even) id,
                        // because a drop removes both the direct and the
                        // RC-induced occurrence for this read.
                        const Shasta2AnchorId canon = v[t].second & ~Shasta2AnchorId(1);
                        journeyTieDrops.push_back({canon, ReadId(readIdValue)});
                        if(remainingMembers[canon] > 0) remainingMembers[canon]--;
                    }
                    deduped.push_back({v[i].first, keeper});
                }
                i = j;
            }
            v.swap(deduped);
        }

        for(const auto& [pos, aid] : v) { static_cast<void>(pos); if(isHetAnchor(aid)) hetOccurrences++; }

        std::vector<Shasta2AnchorId>& out = filteredJourneys[2 * readIdValue];
        std::vector<Shasta2AnchorId>& outRc = filteredJourneys[2 * readIdValue + 1];
        out.reserve(v.size());
        for(const auto& [position, anchorId] : v) {
            static_cast<void>(position);
            out.push_back(anchorId);
        }
        outRc.resize(out.size());
        for(size_t k = 0; k < out.size(); k++) {
            outRc[k] = out[out.size() - 1 - k] ^ 1ULL;
        }
    }
    strand0.clear();
    strand0.shrink_to_fit();
    if(tieGroups > 0) {
        cout << timestamp << "Journeys rebuild: resolved " << tieGroups
             << " same-position tie(s) across independent anchors, dropping "
             << unitsDropped << " occurrence(s) (one survivor kept each)." << endl;
        cout << timestamp << "  of the dropped occurrences, " << hetDropped
             << " were het anchors and " << primaryDropped << " primary; "
             << hetOccurrences << " het occurrences survive in journeys ("
             << (hetOccurrences + hetDropped == 0 ? 0.0 :
                 100.0 * double(hetDropped) / double(hetOccurrences + hetDropped))
             << "% of het occurrences lost to ties)." << endl;
        uint64_t evicted = 0;
        for(uint64_t a = 0; a < anchorCount; a++) {
            if(initialOccurrences[a] > 0 && remainingOccurrences[a] == 0) ++evicted;
        }
        if(hetDropped > 0) {
            cout << timestamp << "  het occurrences lost because the colliding "
                    "anchor was at its coverage floor: " << hetLostToFloor
                 << ", for any other reason: " << hetLostOther << endl;
            if(!floorWinnerOriginalCoverage.empty()) {
                std::sort(floorWinnerOriginalCoverage.begin(),
                          floorWinnerOriginalCoverage.end());
                cout << timestamp << "    original coverage of those winners:";
                for(uint64_t i = 0; i < floorWinnerOriginalCoverage.size(); ) {
                    uint64_t j = i;
                    while(j < floorWinnerOriginalCoverage.size() &&
                          floorWinnerOriginalCoverage[j] == floorWinnerOriginalCoverage[i]) j++;
                    cout << " " << floorWinnerOriginalCoverage[i] << "x:" << (j - i);
                    i = j;
                }
                cout << endl;
                cout << timestamp << "    tie groups where a het lost, by size:";
                std::sort(floorGroupSize.begin(), floorGroupSize.end());
                for(uint64_t i = 0; i < floorGroupSize.size(); ) {
                    uint64_t j = i;
                    while(j < floorGroupSize.size() && floorGroupSize[j] == floorGroupSize[i]) j++;
                    cout << " " << floorGroupSize[i] << "-way:" << (j - i);
                    i = j;
                }
                cout << " (het in group: ";
                std::sort(floorGroupHets.begin(), floorGroupHets.end());
                for(uint64_t i = 0; i < floorGroupHets.size(); ) {
                    uint64_t j = i;
                    while(j < floorGroupHets.size() && floorGroupHets[j] == floorGroupHets[i]) j++;
                    cout << floorGroupHets[i] << ":" << (j - i) << " ";
                    i = j;
                }
                cout << ")" << endl;
                cout << timestamp << "    of the displaced primary anchor's members, "
                        "the share also in the het anchor:";
                for(uint64_t i = 0; i < overlapShared.size(); i++) {
                    cout << " " << overlapShared[i] << "/" << overlapPrimary[i];
                }
                cout << endl;
            }
        }
        cout << timestamp << "  anchors evicted from every journey: " << evicted
             << (evicted == 0 ? "  (none -- an anchor's last occurrence always wins)" : "")
             << endl;
    }

    // Pass C: rebuild the journeys VectorOfVectors in place from
    // filteredJourneys.
    journeys.remove();
    journeys.createNew(largeDataName("Shasta2Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        journeys.incrementCount(oidValue, filteredJourneys[oidValue].size());
    }
    journeys.beginPass2();
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const auto journey = journeys[oidValue];
        const std::vector<Shasta2AnchorId>& filtered = filteredJourneys[oidValue];
        DINARA_ASSERT(journey.size() == filtered.size());
        for(uint64_t i = 0; i < filtered.size(); i++) {
            journey[i] = filtered[i];
        }
    }
    journeys.endPass2(false, true);
    filteredJourneys.clear();
    filteredJourneys.shrink_to_fit();

    // Pass D: reconcile positionInJourney for every anchor's marker infos
    // (reset all to invalid, then set from the rebuilt journeys).
    for(uint64_t anchorId = 0; anchorId < anchors.anchorMarkerInfos.size(); anchorId++) {
        for(Shasta2AnchorMarkerInfo& markerInfo : anchors.anchorMarkerInfos[anchorId]) {
            markerInfo.positionInJourney = invalid<uint32_t>;
        }
    }
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(oidValue));
        const auto journey = journeys[oidValue];
        for(uint64_t position = 0; position < journey.size(); position++) {
            const Shasta2AnchorId anchorId = journey[position];
            span<Shasta2AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
            const auto it = std::lower_bound(markerInfos.begin(), markerInfos.end(), orientedReadId,
                [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                    return info.orientedReadId < oid;
                });
            if(it != markerInfos.end() and it->orientedReadId == orientedReadId) {
                it->positionInJourney = uint32_t(position);
            }
        }
    }

    performanceLog << timestamp << "Journeys rebuild (new anchors) ends." << endl;
}



// Access from binary data.
Shasta2Journeys::Shasta2Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Shasta2Journeys"));
}
