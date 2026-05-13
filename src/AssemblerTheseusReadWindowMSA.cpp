// AssemblerTheseusReadWindowMSA.cpp
//
// Diagnostic prototype: partition physical reads into disjoint Shasta2 anchor
// windows, then (optionally) run Theseus MSAs on consecutive backbone anchor
// pairs using marker-graph MSA recruitment rules.

#include "Assembler.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"
#include "invalid.hpp"

#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Serializes SITE / diagnostic lines from printReadWindowVariationSitesFromMsa when
// backbone-pair MSAs run on multiple threads.
mutex readWindowVariationSiteLogMutex;

constexpr array<uint64_t, 14> histogramUpperBounds = {
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
    numeric_limits<uint64_t>::max()
};

constexpr bool runTheseusMsa = false;
// Backbone-only consecutive Shasta2 anchor pairs per window (recruitment matches TheseusMGMSA).
constexpr bool runAnchorWindowBackbonePairTheseusMsa = true;
constexpr uint32_t minSharedAnchorsForRescue = 4;
constexpr uint64_t minAnchorCoverageForWindowPairMsa = 6;
constexpr uint64_t maxReadsPerWindowAnchorPair = 1024;
constexpr uint64_t oneSidedOffsetRatioNumerator = 11;
constexpr uint64_t oneSidedOffsetRatioDenominator = 10;

// Restrict anchor-window planning to one strand-0 backbone read (diagnostic).
constexpr bool restrictAnchorWindowPlannerToSingleBackboneReadId = true;
constexpr ReadId anchorWindowPlannerOnlyBackboneReadId = ReadId(3729);

// Variation sites: same thresholds as AssemblerTheseusMarkerGraphMSA (ref and alt support >= 3).
constexpr uint64_t rwMinSnpRefSupport = 3;
constexpr uint64_t rwMinSnpAltSupport = 3;
constexpr uint64_t rwMinReportedAltLength = 16;
constexpr uint64_t rwMinFilteredHomopolymerRunLength = 3;

struct ReadWindowTask {
    uint32_t windowId = 0;
    ReadId backboneReadId;
    vector<OrientedReadId> orientedReads;
    vector<ReadId> claimedReads;
    vector<uint32_t> alignmentIds;
};

struct EvidenceOccurrence {
    uint32_t windowId = 0;
    uint32_t row = 0;
    uint32_t alignmentId = 0;
};

struct CandidateEvidence {
    OrientedReadId orientedReadId;
    uint32_t alignmentId = 0;
};

struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t sharedBackboneAnchors = 0;
};

struct AnchorWindowTask {
    uint32_t windowId = 0;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin = 0;
    uint32_t backboneEnd = 0;
    uint64_t claimedAnchorCount = 0;
    vector<AnchorWindowReadInterval> readIntervals;
};

struct AnchorWindowCandidate {
    OrientedReadId backboneOrientedReadId;
    uint32_t begin = 0;
    uint32_t end = 0;
    uint64_t baseSpan = 0;
    uint64_t readLength = 0;
    uint32_t generation = 0;
};

struct AnchorWindowCandidateLess {
    bool operator()(const AnchorWindowCandidate& a, const AnchorWindowCandidate& b) const
    {
        if(a.baseSpan != b.baseSpan) {
            return a.baseSpan < b.baseSpan;
        }
        const uint32_t anchorCountA = a.end - a.begin;
        const uint32_t anchorCountB = b.end - b.begin;
        if(anchorCountA != anchorCountB) {
            return anchorCountA < anchorCountB;
        }
        if(a.readLength != b.readLength) {
            return a.readLength < b.readLength;
        }
        return a.backboneOrientedReadId.getReadId() > b.backboneOrientedReadId.getReadId();
    }
};

struct ThreadCounters {
    uint64_t windows = 0;
    uint64_t skippedSmallWindows = 0;
    uint64_t rows = 0;
    uint64_t bases = 0;
    double msaSeconds = 0.;
};

template<class T> void printWrappedItems(
    ostream& s,
    const string& prefix,
    const vector<T>& items,
    uint64_t itemsPerLine)
{
    s << timestamp << prefix;
    for(uint64_t i=0; i<items.size(); i++) {
        if((i % itemsPerLine) == 0) {
            s << "\n" << timestamp << "  ";
        }
        s << items[i] << " ";
    }
    s << endl;
}

void addToHistogram(array<uint64_t, histogramUpperBounds.size()>& histogram, uint64_t value)
{
    for(uint64_t i=0; i<histogramUpperBounds.size(); i++) {
        if(value <= histogramUpperBounds[i]) {
            ++histogram[i];
            return;
        }
    }
}

string histogramToString(const array<uint64_t, histogramUpperBounds.size()>& histogram)
{
    ostringstream s;
    uint64_t previous = 0;
    for(uint64_t i=0; i<histogram.size(); i++) {
        if(i != 0) {
            s << ",";
        }
        const uint64_t upper = histogramUpperBounds[i];
        if(upper == numeric_limits<uint64_t>::max()) {
            s << ">" << previous;
        } else if(previous + 1 == upper) {
            s << upper;
        } else {
            s << (previous + 1) << "-" << upper;
        }
        s << ":" << histogram[i];
        previous = upper;
    }
    return s.str();
}

string extractWholeOrientedReadSequence(const Reads& reads, OrientedReadId oid)
{
    const uint32_t length = uint32_t(reads.getRead(oid.getReadId()).baseCount);
    string sequence;
    sequence.reserve(length);
    for(uint32_t pos=0; pos<length; pos++) {
        sequence.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return sequence;
}

struct MsaSegment {
    string sequence;
    uint32_t begin = 0;
    uint32_t end = 0;
};

struct MsaSequenceInfo {
    OrientedReadId oid;
    string sequence;
    uint32_t begin = 0;
    uint32_t end = 0;
    bool hasBothAnchors = false;
    bool isEndsFree = false;
    char anchorSide = 'B';
};

MsaSegment extractMsaSegmentFromOrdinals(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    MsaSegment segment;
    if(ordinalA == ordinalB) {
        return segment;
    }

    const uint32_t ord0 = min(ordinalA, ordinalB);
    const uint32_t ord1 = max(ordinalA, ordinalB);
    const auto readMarkers = markers[oid.getValue()];
    if(ord1 >= readMarkers.size()) {
        return segment;
    }

    const uint32_t kHalf = uint32_t(k / 2);
    segment.begin = readMarkers[ord0].position + kHalf;
    segment.end = readMarkers[ord1].position + kHalf;
    if(segment.end <= segment.begin) {
        return MsaSegment{};
    }

    segment.sequence.reserve(segment.end - segment.begin);
    for(uint32_t pos=segment.begin; pos<segment.end; pos++) {
        segment.sequence.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return segment;
}

MsaSegment extractMsaSegmentFromBases(
    const Reads& reads,
    OrientedReadId oid,
    uint32_t begin,
    uint32_t end)
{
    MsaSegment segment;
    const uint32_t readLength = uint32_t(reads.getRead(oid.getReadId()).baseCount);
    begin = min(begin, readLength);
    end = min(end, readLength);
    if(end <= begin) {
        return segment;
    }

    segment.begin = begin;
    segment.end = end;
    segment.sequence.reserve(end - begin);
    for(uint32_t pos=begin; pos<end; pos++) {
        segment.sequence.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return segment;
}

bool msaSegmentContainsMarkerAtLeft(
    const MsaSegment& segment,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinal)
{
    const auto readMarkers = markers[oid.getValue()];
    if(ordinal >= readMarkers.size()) {
        return false;
    }
    const uint32_t anchorPosition = readMarkers[ordinal].position + uint32_t(k / 2);
    return segment.begin == anchorPosition && segment.end > anchorPosition;
}

bool msaSegmentContainsMarkerAtRight(
    const MsaSegment& segment,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinal)
{
    const auto readMarkers = markers[oid.getValue()];
    if(ordinal >= readMarkers.size()) {
        return false;
    }
    const uint32_t anchorPosition = readMarkers[ordinal].position + uint32_t(k / 2);
    return segment.begin < anchorPosition && segment.end == anchorPosition;
}

struct BackbonePairJob {
    uint32_t windowIndex = 0;
    uint32_t backboneJ = 0;
};

// Built on worker threads (recruitment + sequence extraction). Theseus runs serially.
struct PreparedBackbonePairMsaJob {
    bool ready = false;
    BackbonePairJob job;
    Shasta2AnchorId leftId = 0;
    Shasta2AnchorId rightId = 0;
    vector<MsaSequenceInfo> sequenceInfos;
    uint32_t focalWindowBegin = 0;
    uint32_t focalWindowEnd = 0;
    uint32_t focalPaddedBegin = 0;
    uint32_t focalPaddedEnd = 0;
};

struct BackbonePairMsaCounters {
    uint64_t pairJobsScheduled = 0;
    uint64_t skippedCoverage = 0;
    uint64_t skippedOrdinal = 0;
    uint64_t skippedEmptyFocalWindow = 0;
    uint64_t skippedSingleSequence = 0;
    uint64_t skippedReverseOrderBothAnchorReads = 0;
    uint64_t msasRun = 0;
    uint64_t leftOnlySequences = 0;
    uint64_t rightOnlySequences = 0;
    uint64_t failedLeftOnlyAnchorChecks = 0;
    uint64_t failedRightOnlyAnchorChecks = 0;
    double msaCpuSeconds = 0.;
};

void mergeBackbonePairMsaCounters(BackbonePairMsaCounters& total, const BackbonePairMsaCounters& part)
{
    total.skippedCoverage += part.skippedCoverage;
    total.skippedOrdinal += part.skippedOrdinal;
    total.skippedEmptyFocalWindow += part.skippedEmptyFocalWindow;
    total.skippedSingleSequence += part.skippedSingleSequence;
    total.skippedReverseOrderBothAnchorReads += part.skippedReverseOrderBothAnchorReads;
    total.msasRun += part.msasRun;
    total.leftOnlySequences += part.leftOnlySequences;
    total.rightOnlySequences += part.rightOnlySequences;
    total.failedLeftOnlyAnchorChecks += part.failedLeftOnlyAnchorChecks;
    total.failedRightOnlyAnchorChecks += part.failedRightOnlyAnchorChecks;
    total.msaCpuSeconds += part.msaCpuSeconds;
}

string rwDisplayAllele(const string& allele)
{
    constexpr size_t maxShown = 24;
    if(allele.empty()) {
        return "-";
    }
    if(allele.size() <= maxShown) {
        return allele;
    }
    return allele.substr(0, maxShown) + "...(" + to_string(allele.size()) + "bp)";
}

bool rwIsCanonicalBase(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

vector<string> rwParseMsaFasta(const string& text)
{
    vector<string> sequences;
    string current;
    istringstream in(text);
    string line;
    while(getline(in, line)) {
        if(line.empty()) {
            continue;
        }
        if(line[0] == '>') {
            if(!current.empty()) {
                sequences.push_back(current);
                current.clear();
            }
        } else {
            current += line;
        }
    }
    if(!current.empty()) {
        sequences.push_back(current);
    }
    return sequences;
}

string rwOrientedReadList(const vector<uint64_t>& indexes, const vector<MsaSequenceInfo>& sequenceInfos)
{
    string s;
    for(size_t i=0; i<indexes.size(); i++) {
        if(i) {
            s += ",";
        }
        s += sequenceInfos[indexes[i]].oid.getString();
    }
    return s;
}

void rwAddReadIndex(vector<uint64_t>& indexes, uint64_t index)
{
    if(find(indexes.begin(), indexes.end(), index) == indexes.end()) {
        indexes.push_back(index);
    }
}

string rwAlleleType(const string& ref, const string& alt)
{
    if(ref.size() == 1 && alt.size() == 1) {
        return "SNP";
    }
    if(ref.size() == alt.size()) {
        return "MNP";
    }
    if(ref.size() < alt.size()) {
        return "INS";
    }
    if(ref.size() > alt.size()) {
        return "DEL";
    }
    return "COMPLEX";
}

bool rwShouldReportAllele(const string& type, const string& alt)
{
    return type == "SNP" || alt.size() >= rwMinReportedAltLength;
}

bool rwIsHomopolymerAt(const string& seq, uint32_t pos)
{
    if(pos >= seq.size()) {
        return false;
    }

    const size_t p = size_t(pos);
    const char base = seq[p];
    if(!rwIsCanonicalBase(base)) {
        return false;
    }

    size_t begin = p;
    while(begin > 0 && seq[begin - 1] == base) {
        --begin;
    }
    size_t end = p + 1;
    while(end < seq.size() && seq[end] == base) {
        ++end;
    }
    return end - begin >= rwMinFilteredHomopolymerRunLength;
}

bool rwSnpTouchesHomopolymer(const string& focalTargetSequence, uint32_t pos, char altBase)
{
    if(rwIsHomopolymerAt(focalTargetSequence, pos)) {
        return true;
    }
    if(pos >= focalTargetSequence.size() || !rwIsCanonicalBase(altBase)) {
        return false;
    }
    string altTargetSequence = focalTargetSequence;
    altTargetSequence[pos] = altBase;
    return rwIsHomopolymerAt(altTargetSequence, pos);
}

void printReadWindowVariationSitesFromMsa(
    uint64_t pairLabel,
    Shasta2AnchorId leftAnchorId,
    Shasta2AnchorId rightAnchorId,
    const vector<MsaSequenceInfo>& sequenceInfos,
    const vector<string>& alignedSequences,
    uint32_t focalBegin,
    uint32_t focalEnd,
    uint32_t reportBegin,
    uint32_t reportEnd,
    double msaSeconds)
{
    if(focalEnd <= focalBegin) {
        return;
    }

    if(alignedSequences.size() != sequenceInfos.size() || alignedSequences.empty()) {
        lock_guard<mutex> lock(readWindowVariationSiteLogMutex);
        cout << timestamp << "[TheseusReadWindowMSA] pairLabel=" << pairLabel
             << " skipped: MSA sequence count mismatch." << endl;
        return;
    }

    const size_t columnCount = alignedSequences.front().size();
    for(const string& s: alignedSequences) {
        if(s.size() != columnCount) {
            lock_guard<mutex> lock(readWindowVariationSiteLogMutex);
            cout << timestamp << "[TheseusReadWindowMSA] pairLabel=" << pairLabel
                 << " skipped: MSA rows have inconsistent lengths." << endl;
            return;
        }
    }

    if(columnCount == 0) {
        return;
    }

    vector<size_t> firstNonGap(alignedSequences.size(), columnCount);
    vector<size_t> lastNonGap(alignedSequences.size(), 0);
    vector<uint8_t> hasBase(alignedSequences.size(), 0);
    for(size_t r=0; r<alignedSequences.size(); r++) {
        for(size_t c=0; c<columnCount; c++) {
            if(alignedSequences[r][c] != '-') {
                firstNonGap[r] = min(firstNonGap[r], c);
                lastNonGap[r] = max(lastNonGap[r], c);
                hasBase[r] = 1;
            }
        }
    }

    vector<uint32_t> targetPosByColumn(columnCount, invalid<uint32_t>);
    vector<size_t> columnByTargetOffset;
    columnByTargetOffset.reserve(focalEnd - focalBegin);
    uint32_t focalPos = focalBegin;
    for(size_t c=0; c<columnCount; c++) {
        const char refBase = alignedSequences[0][c];
        if(refBase == '-') {
            continue;
        }
        targetPosByColumn[c] = focalPos++;
        columnByTargetOffset.push_back(c);
    }

    vector<size_t> targetOffsetByColumn(columnCount, invalid<size_t>);
    for(size_t offset=0; offset<columnByTargetOffset.size(); offset++) {
        targetOffsetByColumn[columnByTargetOffset[offset]] = offset;
    }
    string focalTargetSequence;
    focalTargetSequence.reserve(columnByTargetOffset.size());
    for(const size_t c: columnByTargetOffset) {
        const char base = alignedSequences[0][c];
        focalTargetSequence.push_back(rwIsCanonicalBase(base) ? base : 'N');
    }

    vector<uint8_t> isDirtyTargetOffset(columnByTargetOffset.size(), 0);
    for(size_t offset=0; offset<columnByTargetOffset.size(); offset++) {
        const size_t c = columnByTargetOffset[offset];
        const char refBase = alignedSequences[0][c];
        if(!rwIsCanonicalBase(refBase)) {
            continue;
        }
        for(size_t r=1; r<alignedSequences.size(); r++) {
            if(!hasBase[r] || c < firstNonGap[r] || c > lastNonGap[r]) {
                continue;
            }
            const char base = alignedSequences[r][c];
            if(base == '-') {
                isDirtyTargetOffset[offset] = 1;
            } else if(rwIsCanonicalBase(base) && base != refBase) {
                isDirtyTargetOffset[offset] = 1;
            }
        }
    }

    size_t gapRunBegin = invalid<size_t>;
    for(size_t c=0; c<=columnCount; c++) {
        const bool isFocalGap = (c < columnCount && alignedSequences[0][c] == '-');
        if(isFocalGap && gapRunBegin == invalid<size_t>) {
            gapRunBegin = c;
        }
        if((!isFocalGap || c == columnCount) && gapRunBegin != invalid<size_t>) {
            const size_t gapRunEnd = c;
            size_t anchorColumn = invalid<size_t>;
            for(size_t j=gapRunBegin; j>0; --j) {
                if(alignedSequences[0][j - 1] != '-') {
                    anchorColumn = j - 1;
                    break;
                }
            }
            if(anchorColumn != invalid<size_t> &&
               targetPosByColumn[anchorColumn] != invalid<uint32_t> &&
               targetOffsetByColumn[anchorColumn] != invalid<size_t>) {
                for(size_t r=1; r<alignedSequences.size(); r++) {
                    if(!hasBase[r] || anchorColumn < firstNonGap[r] || anchorColumn > lastNonGap[r]) {
                        continue;
                    }
                    for(size_t j=gapRunBegin; j<gapRunEnd; j++) {
                        if(rwIsCanonicalBase(alignedSequences[r][j])) {
                            isDirtyTargetOffset[targetOffsetByColumn[anchorColumn]] = 1;
                            break;
                        }
                    }
                }
            }
            gapRunBegin = invalid<size_t>;
        }
    }

    using AlleleKey = tuple<string, string, string>;
    for(size_t beginOffset=0; beginOffset<isDirtyTargetOffset.size();) {
        if(!isDirtyTargetOffset[beginOffset]) {
            ++beginOffset;
            continue;
        }
        size_t endOffset = beginOffset;
        while(endOffset + 1 < isDirtyTargetOffset.size() &&
              isDirtyTargetOffset[endOffset + 1]) {
            ++endOffset;
        }
        const uint32_t targetPos = focalBegin + uint32_t(beginOffset);
        const uint32_t targetEnd = focalBegin + uint32_t(endOffset + 1);
        if(targetPos < reportBegin || targetEnd > reportEnd) {
            beginOffset = endOffset + 1;
            continue;
        }

        const size_t beginColumn = columnByTargetOffset[beginOffset];
        size_t endColumn = columnByTargetOffset[endOffset];
        for(size_t c=endColumn + 1; c<columnCount && alignedSequences[0][c] == '-'; c++) {
            endColumn = c;
        }

        string ref;
        for(size_t c=beginColumn; c<=endColumn; c++) {
            const char base = alignedSequences[0][c];
            if(rwIsCanonicalBase(base)) {
                ref.push_back(base);
            }
        }
        if(ref.empty()) {
            beginOffset = endOffset + 1;
            continue;
        }

        vector<uint64_t> refReads;
        map<AlleleKey, vector<uint64_t>> altReadsByAllele;
        for(size_t r=0; r<alignedSequences.size(); r++) {
            if(!hasBase[r] || beginColumn < firstNonGap[r] || endColumn > lastNonGap[r]) {
                continue;
            }
            string allele;
            for(size_t c=beginColumn; c<=endColumn; c++) {
                const char base = alignedSequences[r][c];
                if(rwIsCanonicalBase(base)) {
                    allele.push_back(base);
                }
            }
            if(allele == ref) {
                refReads.push_back(r);
            } else {
                const AlleleKey key{rwAlleleType(ref, allele), ref, allele};
                rwAddReadIndex(altReadsByAllele[key], r);
            }
        }

        uint64_t maxAltSupport = 0;
        map<AlleleKey, vector<uint64_t>> reportableAltReadsByAllele;
        for(const auto& [key, altReads]: altReadsByAllele) {
            const auto& [type, refA, alt] = key;
            (void) refA;
            if(!rwShouldReportAllele(type, alt)) {
                continue;
            }
            if(type == "SNP" &&
               rwSnpTouchesHomopolymer(focalTargetSequence, uint32_t(beginOffset), alt[0])) {
                continue;
            }
            reportableAltReadsByAllele[key] = altReads;
            maxAltSupport = max<uint64_t>(maxAltSupport, altReads.size());
        }
        if(refReads.size() < rwMinSnpRefSupport || maxAltSupport < rwMinSnpAltSupport) {
            beginOffset = endOffset + 1;
            continue;
        }

        {
            lock_guard<mutex> lock(readWindowVariationSiteLogMutex);
            cout << timestamp << "[TheseusReadWindowMSA] SITE"
                 << " pairLabel=" << pairLabel
                 << " anchors=" << leftAnchorId << "->" << rightAnchorId
                 << " focal=" << sequenceInfos.front().oid
                 << " pos=" << targetPos
                 << " ref=" << rwDisplayAllele(ref)
                 << " refN=" << refReads.size()
                 << " refReads=" << rwOrientedReadList(refReads, sequenceInfos)
                 << " msaSeconds=" << fixed << setprecision(6) << msaSeconds;
            for(const auto& [key, altReads]: reportableAltReadsByAllele) {
                const auto& [type, refK, alt] = key;
                (void) refK;
                cout << " alt=" << type << ":" << rwDisplayAllele(ref) << ">"
                     << rwDisplayAllele(alt)
                     << ":N=" << altReads.size()
                     << ":reads=" << rwOrientedReadList(altReads, sequenceInfos);
            }
            cout << defaultfloat << endl;
        }

        beginOffset = endOffset + 1;
    }
}

void prepareBackbonePairMsaJob(
    const Assembler& assembler,
    const vector<AnchorWindowTask>& anchorWindows,
    const Shasta2Anchors& shasta2Anchors,
    const Shasta2Journeys& shasta2Journeys,
    const BackbonePairJob& job,
    uint64_t orientedReadCount,
    PreparedBackbonePairMsaJob& prep,
    BackbonePairMsaCounters& c)
{
    prep = PreparedBackbonePairMsaJob{};
    DINARA_ASSERT(assembler.markers);
    const Reads& reads = assembler.getReads();
    const auto& markers = *assembler.markers;
    const uint64_t k = assembler.assemblerInfo->k;

    const AnchorWindowTask& task = anchorWindows[job.windowIndex];
    const OrientedReadId focalOid = task.backboneOrientedReadId;
    const auto journey = shasta2Journeys[focalOid];
    const uint32_t j = job.backboneJ;
    if(j + 1 >= journey.size()) {
        ++c.skippedOrdinal;
        return;
    }

    const Shasta2AnchorId leftId = journey[j];
    const Shasta2AnchorId rightId = journey[j + 1];
    const Shasta2Anchor leftAnchor = shasta2Anchors[leftId];
    const Shasta2Anchor rightAnchor = shasta2Anchors[rightId];
    if(leftAnchor.size() < minAnchorCoverageForWindowPairMsa ||
        rightAnchor.size() < minAnchorCoverageForWindowPairMsa) {
        ++c.skippedCoverage;
        return;
    }

    const uint32_t focalLeftOrdinal = shasta2Anchors.getOrdinal(leftId, focalOid);
    const uint32_t focalRightOrdinal = shasta2Anchors.getOrdinal(rightId, focalOid);
    if(focalLeftOrdinal == invalid<uint32_t> || focalRightOrdinal == invalid<uint32_t>) {
        ++c.skippedOrdinal;
        return;
    }

    MsaSegment focalWindowSegment = extractMsaSegmentFromOrdinals(
        reads, markers, k, focalOid, focalLeftOrdinal, focalRightOrdinal);
    if(focalWindowSegment.sequence.empty()) {
        ++c.skippedEmptyFocalWindow;
        return;
    }

    unordered_map<uint64_t, uint32_t> leftOrdinals;
    unordered_map<uint64_t, uint32_t> rightOrdinals;
    for(const Shasta2AnchorMarkerInfo& info: leftAnchor) {
        leftOrdinals.try_emplace(info.orientedReadId.getValue(), info.ordinal);
    }
    for(const Shasta2AnchorMarkerInfo& info: rightAnchor) {
        rightOrdinals.try_emplace(info.orientedReadId.getValue(), info.ordinal);
    }

    vector<uint64_t> candidateValues;
    candidateValues.reserve(leftOrdinals.size() + rightOrdinals.size());
    for(const auto& [oidValue, ordinal]: leftOrdinals) {
        (void) ordinal;
        if(oidValue != focalOid.getValue()) {
            candidateValues.push_back(oidValue);
        }
    }
    for(const auto& [oidValue, ordinal]: rightOrdinals) {
        (void) ordinal;
        if(oidValue != focalOid.getValue() && !leftOrdinals.contains(oidValue)) {
            candidateValues.push_back(oidValue);
        }
    }
    sort(candidateValues.begin(), candidateValues.end());
    stable_sort(candidateValues.begin(), candidateValues.end(),
        [&](uint64_t oidValueA, uint64_t oidValueB) {
            const bool aHasBoth =
                leftOrdinals.contains(oidValueA) && rightOrdinals.contains(oidValueA);
            const bool bHasBoth =
                leftOrdinals.contains(oidValueB) && rightOrdinals.contains(oidValueB);
            if(aHasBoth != bHasBoth) {
                return aHasBoth;
            }
            return oidValueA < oidValueB;
        });

    uint64_t approximateSpanSum = 0;
    uint64_t approximateSpanCount = 0;
    for(const auto& [oidValue, leftOrdinal]: leftOrdinals) {
        auto itRight = rightOrdinals.find(oidValue);
        if(itRight == rightOrdinals.end()) {
            continue;
        }
        if(oidValue >= orientedReadCount) {
            continue;
        }
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        const auto readMarkers = markers[oid.getValue()];
        const uint32_t rightOrdinal = itRight->second;
        if(leftOrdinal >= readMarkers.size() || rightOrdinal >= readMarkers.size()) {
            continue;
        }
        if(leftOrdinal >= rightOrdinal) {
            continue;
        }
        const uint32_t kHalf = uint32_t(k / 2);
        const uint32_t begin = readMarkers[leftOrdinal].position + kHalf;
        const uint32_t end = readMarkers[rightOrdinal].position + kHalf;
        if(end > begin) {
            approximateSpanSum += (end - begin);
            ++approximateSpanCount;
        }
    }
    if(approximateSpanCount == 0) {
        approximateSpanSum = focalWindowSegment.end - focalWindowSegment.begin;
        approximateSpanCount = 1;
    }
    const uint64_t averageSpan = max<uint64_t>(1, approximateSpanSum / approximateSpanCount);
    const uint32_t approximateSpan = uint32_t(max<uint64_t>(
        1,
        (averageSpan * oneSidedOffsetRatioNumerator +
            oneSidedOffsetRatioDenominator - 1) /
            oneSidedOffsetRatioDenominator));
    const uint32_t focalWindowLength = focalWindowSegment.end - focalWindowSegment.begin;
    const uint32_t alignmentPadding =
        (approximateSpan > focalWindowLength) ? (approximateSpan - focalWindowLength) : 0;
    MsaSegment focalSegment = extractMsaSegmentFromBases(
        reads,
        focalOid,
        (focalWindowSegment.begin > alignmentPadding) ?
            (focalWindowSegment.begin - alignmentPadding) : 0,
        focalWindowSegment.end + alignmentPadding);
    if(focalSegment.sequence.empty()) {
        ++c.skippedEmptyFocalWindow;
        return;
    }

    vector<MsaSequenceInfo> sequenceInfos;
    sequenceInfos.reserve(size_t(maxReadsPerWindowAnchorPair));
    sequenceInfos.push_back(MsaSequenceInfo{
        focalOid,
        focalSegment.sequence,
        focalSegment.begin,
        focalSegment.end,
        true,
        false,
        'B'});

    for(const uint64_t oidValue64: candidateValues) {
        if(sequenceInfos.size() >= maxReadsPerWindowAnchorPair) {
            break;
        }
        if(oidValue64 >= orientedReadCount) {
            continue;
        }

        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue64));
        const bool hasLeft = leftOrdinals.contains(oidValue64);
        const bool hasRight = rightOrdinals.contains(oidValue64);
        MsaSegment segment;
        bool hasBothAnchors = false;
        bool isEndsFree = true;
        char anchorSide = 'B';

        if(hasLeft && hasRight) {
            if(leftOrdinals[oidValue64] >= rightOrdinals[oidValue64]) {
                ++c.skippedReverseOrderBothAnchorReads;
                continue;
            }
            const auto readMarkers = markers[oid.getValue()];
            const uint32_t leftOrdinal = leftOrdinals[oidValue64];
            const uint32_t rightOrdinal = rightOrdinals[oidValue64];
            if(leftOrdinal < readMarkers.size() && rightOrdinal < readMarkers.size()) {
                const uint32_t kHalf = uint32_t(k / 2);
                const uint32_t windowBegin = readMarkers[leftOrdinal].position + kHalf;
                const uint32_t windowEnd = readMarkers[rightOrdinal].position + kHalf;
                if(windowEnd > windowBegin) {
                    segment = extractMsaSegmentFromBases(
                        reads,
                        oid,
                        (windowBegin > alignmentPadding) ? (windowBegin - alignmentPadding) : 0,
                        windowEnd + alignmentPadding);
                }
            }
            hasBothAnchors = true;
            isEndsFree = false;
            anchorSide = 'B';
        }

        if(segment.sequence.empty() && hasLeft) {
            const auto readMarkers = markers[oid.getValue()];
            const uint32_t leftOrdinal = leftOrdinals[oidValue64];
            if(leftOrdinal < readMarkers.size()) {
                const uint32_t begin =
                    readMarkers[leftOrdinal].position + uint32_t(k / 2);
                segment = extractMsaSegmentFromBases(
                    reads, oid, begin, begin + approximateSpan);
            }
            anchorSide = 'L';
        }

        if(segment.sequence.empty() && hasRight) {
            const uint32_t rightOrdinal = rightOrdinals[oidValue64];
            const auto readMarkers = markers[oid.getValue()];
            if(rightOrdinal < readMarkers.size()) {
                const uint32_t end =
                    readMarkers[rightOrdinal].position + uint32_t(k / 2);
                const uint32_t begin = (end > approximateSpan) ? (end - approximateSpan) : 0;
                segment = extractMsaSegmentFromBases(reads, oid, begin, end);
            }
            anchorSide = 'R';
        }

        if(segment.sequence.empty()) {
            continue;
        }

        if(anchorSide == 'L') {
            if(!msaSegmentContainsMarkerAtLeft(
                segment, markers, k, oid, leftOrdinals[oidValue64])) {
                ++c.failedLeftOnlyAnchorChecks;
                continue;
            }
            ++c.leftOnlySequences;
        } else if(anchorSide == 'R') {
            if(!msaSegmentContainsMarkerAtRight(
                segment, markers, k, oid, rightOrdinals[oidValue64])) {
                ++c.failedRightOnlyAnchorChecks;
                continue;
            }
            ++c.rightOnlySequences;
        }

        sequenceInfos.push_back(MsaSequenceInfo{
            oid,
            segment.sequence,
            segment.begin,
            segment.end,
            hasBothAnchors,
            isEndsFree,
            anchorSide});
    }

    if(sequenceInfos.size() < 2) {
        ++c.skippedSingleSequence;
        return;
    }

    prep.ready = true;
    prep.job = job;
    prep.leftId = leftId;
    prep.rightId = rightId;
    prep.sequenceInfos = std::move(sequenceInfos);
    prep.focalWindowBegin = focalWindowSegment.begin;
    prep.focalWindowEnd = focalWindowSegment.end;
    prep.focalPaddedBegin = focalSegment.begin;
    prep.focalPaddedEnd = focalSegment.end;
}

void runTheseusBackbonePairOnPrepared(
    const PreparedBackbonePairMsaJob& prep,
    const vector<AnchorWindowTask>& anchorWindows,
    BackbonePairMsaCounters& c)
{
    if(!prep.ready) {
        return;
    }
    const AnchorWindowTask& task = anchorWindows[prep.job.windowIndex];
    const vector<MsaSequenceInfo>& sequenceInfos = prep.sequenceInfos;

    const auto msaBegin = chrono::steady_clock::now();
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics(false, false);
    theseus::TheseusMSA aligner(
        penalties,
        heuristics,
        sequenceInfos.front().sequence,
        1,
        false);

    for(size_t i=1; i<sequenceInfos.size(); i++) {
        aligner.align(sequenceInfos[i].sequence, 1, false, sequenceInfos[i].isEndsFree);
    }

    ostringstream msaOut;
    aligner.print_as_msa(msaOut);

    const auto msaEnd = chrono::steady_clock::now();
    const double msaSecondsOne = chrono::duration<double>(msaEnd - msaBegin).count();
    c.msaCpuSeconds += msaSecondsOne;
    ++c.msasRun;

    const string msaText = msaOut.str();
    const vector<string> alignedSequences = rwParseMsaFasta(msaText);
    const uint64_t pairLabel =
        (uint64_t(task.windowId) << 32) | uint64_t(prep.job.backboneJ);
    printReadWindowVariationSitesFromMsa(
        pairLabel,
        prep.leftId,
        prep.rightId,
        sequenceInfos,
        alignedSequences,
        prep.focalPaddedBegin,
        prep.focalPaddedEnd,
        prep.focalWindowBegin,
        prep.focalWindowEnd,
        msaSecondsOne);
}

void runShasta2BackbonePairTheseusMsas(
    const Assembler& assembler,
    const vector<AnchorWindowTask>& anchorWindows,
    const Shasta2Anchors& shasta2Anchors,
    const Shasta2Journeys& shasta2Journeys,
    uint64_t threadCount,
    uint64_t orientedReadCount,
    BackbonePairMsaCounters& total)
{
    vector<BackbonePairJob> jobs;
    for(uint32_t wi=0; wi<uint32_t(anchorWindows.size()); wi++) {
        const AnchorWindowTask& task = anchorWindows[wi];
        if(task.backboneBegin + 1 >= task.backboneEnd) {
            continue;
        }
        const auto journey = shasta2Journeys[task.backboneOrientedReadId];
        if(journey.size() < 2) {
            continue;
        }
        const uint32_t pairEndExclusive =
            min(task.backboneEnd - 1, uint32_t(journey.size() - 1));
        for(uint32_t j=task.backboneBegin; j<pairEndExclusive; j++) {
            jobs.push_back(BackbonePairJob{wi, j});
        }
    }

    total.pairJobsScheduled = jobs.size();
    if(jobs.empty()) {
        return;
    }

    // Each thread owns a disjoint subset of consecutive-pair jobs (strided by thread id)
    // and runs prepare + Theseus for every job in that subset (full job per thread).
    threadCount = max<uint64_t>(1, threadCount);
    const uint64_t workerCount = min<uint64_t>(threadCount, jobs.size());
    vector<BackbonePairMsaCounters> threadTotals(workerCount);
    vector<thread> workers;
    workers.reserve(workerCount);

    for(uint64_t t=0; t<workerCount; t++) {
        workers.emplace_back([&, t]() {
            BackbonePairMsaCounters local;
            for(uint64_t i=t; i<jobs.size(); i+=workerCount) {
                PreparedBackbonePairMsaJob prep;
                prepareBackbonePairMsaJob(
                    assembler,
                    anchorWindows,
                    shasta2Anchors,
                    shasta2Journeys,
                    jobs[i],
                    orientedReadCount,
                    prep,
                    local);
                if(prep.ready) {
                    runTheseusBackbonePairOnPrepared(prep, anchorWindows, local);
                }
            }
            threadTotals[t] = std::move(local);
        });
    }
    for(thread& th: workers) {
        th.join();
    }

    for(uint64_t t=0; t<workerCount; t++) {
        mergeBackbonePairMsaCounters(total, threadTotals[t]);
    }
}

} // namespace



void Assembler::computeTheseusReadWindowMSAPrototype(uint64_t threadCount)
{
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    DINARA_ASSERT(assemblerInfo.isOpen);
    DINARA_ASSERT((assemblerInfo->k % 2) == 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    const MappedMemoryOwner shasta2Owner = shasta2MappedMemoryOwner();
    auto shasta2AnchorsLocal = make_shared<Shasta2Anchors>(
        shasta2Owner,
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        threadCount,
        2,
        numeric_limits<uint64_t>::max());
    auto shasta2JourneysLocal = make_shared<Shasta2Journeys>(
        2 * getReads().readCount(),
        shasta2AnchorsLocal,
        threadCount,
        shasta2Owner);
    computeTheseusReadWindowMSAPrototype(shasta2AnchorsLocal, shasta2JourneysLocal, threadCount);
}



void Assembler::computeTheseusReadWindowMSAPrototype(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    uint64_t threadCount)
{
    cout << timestamp << "[TheseusReadWindowMSA] Prototype begins." << endl;
    const auto totalBegin = chrono::steady_clock::now();

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    checkReadGraphIsOpen();
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    DINARA_ASSERT(assemblerInfo.isOpen);
    DINARA_ASSERT((assemblerInfo->k % 2) == 0);
    DINARA_ASSERT(reads->readCount() > 0);
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const uint32_t invalidAlignmentId = numeric_limits<uint32_t>::max();

    {
    // Anchor-interval window prototype.
    // Claim anchor ids, not whole reads. Each accepted window is seeded by a
    // contiguous unclaimed interval on a strand-0 backbone journey, then expanded
    // to contiguous unclaimed intervals on reads touching the seed anchors.
    constexpr uint32_t minBackboneWindowAnchors = 2;
    const uint32_t anchorUnclaimed = numeric_limits<uint32_t>::max();
    const uint64_t anchorCount = shasta2Anchors->size();
    vector<uint32_t> anchorOwner(anchorCount, anchorUnclaimed);
    vector<AnchorWindowTask> anchorWindows;
    vector<uint32_t> touchedEpoch(orientedReadCount, 0);
    vector<uint32_t> touchedMin(orientedReadCount, numeric_limits<uint32_t>::max());
    vector<uint32_t> touchedMax(orientedReadCount, 0);
    vector<uint32_t> touchedCount(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint32_t epoch = 0;
    priority_queue<
        AnchorWindowCandidate,
        vector<AnchorWindowCandidate>,
        AnchorWindowCandidateLess> candidateHeap;
    vector<uint32_t> candidateGeneration(readCount, 0);

    vector<ReadId> anchorReadsByLength;
    anchorReadsByLength.reserve(readCount);
    for(uint64_t readId=0; readId<readCount; readId++) {
        anchorReadsByLength.push_back(ReadId(readId));
    }
    sort(anchorReadsByLength.begin(), anchorReadsByLength.end(),
        [&](ReadId a, ReadId b) {
            const uint64_t lengthA = reads->getRead(a).baseCount;
            const uint64_t lengthB = reads->getRead(b).baseCount;
            if(lengthA != lengthB) {
                return lengthA > lengthB;
            }
            return a < b;
        });

    uint64_t anchorWindowClaimedAnchors = 0;
    uint64_t backboneClaimedAnchors = 0;
    uint64_t nonBackboneClaimedAnchors = 0;
    uint64_t anchorWindowReadIntervals = 0;
    uint64_t anchorWindowSkippedNoJourney = 0;
    uint64_t anchorWindowSkippedShortRuns = 0;
    uint64_t anchorWindowBackboneIntervals = 0;
    uint64_t anchorWindowTouchedReads = 0;
    uint64_t anchorWindowSplitIntervals = 0;
    uint64_t candidateIntervalsPushed = 0;
    uint64_t candidateIntervalsPopped = 0;
    uint64_t staleCandidateIntervals = 0;
    uint64_t discardedOldGenerationCandidates = 0;
    uint64_t maxAnchorWindowReadIntervals = 0;
    uint64_t maxAnchorWindowClaimedAnchors = 0;
    const auto anchorWindowBegin = chrono::steady_clock::now();

    auto intervalBaseSpan = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(begin >= end || oid.getValue() >= markers->size()) {
            return uint64_t(end - begin);
        }
        const auto orientedReadMarkers = (*markers)[oid.getValue()];
        if(orientedReadMarkers.empty()) {
            return uint64_t(end - begin);
        }

        const Shasta2AnchorId leftAnchorId = journey[begin];
        const Shasta2AnchorId rightAnchorId = journey[end - 1];
        const uint32_t leftOrdinal = shasta2Anchors->getOrdinal(leftAnchorId, oid);
        const uint32_t rightOrdinal = shasta2Anchors->getOrdinal(rightAnchorId, oid);
        if(leftOrdinal == invalid<uint32_t> || rightOrdinal == invalid<uint32_t>) {
            return uint64_t(end - begin);
        }

        const uint64_t leftOrdClamped =
            min<uint64_t>(leftOrdinal, orientedReadMarkers.size() - 1);
        const uint64_t leftPosition = orientedReadMarkers[leftOrdClamped].position;
        uint64_t rightPosition = reads->getRead(oid.getReadId()).baseCount;
        if(rightOrdinal < orientedReadMarkers.size()) {
            rightPosition = orientedReadMarkers[rightOrdinal].position;
        }
        return rightPosition > leftPosition ? rightPosition - leftPosition : uint64_t(end - begin);
    };

    auto pushCandidate = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(end - begin < minBackboneWindowAnchors) {
            ++anchorWindowSkippedShortRuns;
            return;
        }
        candidateHeap.push(AnchorWindowCandidate{
            oid,
            begin,
            end,
            intervalBaseSpan(oid, journey, begin, end),
            reads->getRead(oid.getReadId()).baseCount,
            candidateGeneration[uint64_t(oid.getReadId())]});
        ++candidateIntervalsPushed;
    };

    auto pushCurrentUnclaimedIntervals = [&](OrientedReadId oid) {
        if(oid.getValue() >= shasta2Journeys->size()) {
            ++anchorWindowSkippedNoJourney;
            return;
        }
        const auto journey = (*shasta2Journeys)[oid];
        if(journey.empty()) {
            ++anchorWindowSkippedNoJourney;
            return;
        }

        uint32_t position = 0;
        while(position < journey.size()) {
            while(position < journey.size() &&
                  anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                ++position;
            }
            const uint32_t runBegin = position;
            while(position < journey.size() &&
                  anchorOwner[uint64_t(journey[position])] == anchorUnclaimed) {
                ++position;
            }
            const uint32_t runEnd = position;
            if(runBegin != runEnd) {
                pushCandidate(oid, journey, runBegin, runEnd);
            }
        }
    };

    auto createAnchorWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindowTask task;
        task.windowId = windowId;
        task.backboneOrientedReadId = backboneOid;
        task.backboneBegin = seedBegin;
        task.backboneEnd = seedEnd;
        task.readIntervals.push_back(AnchorWindowReadInterval{
            backboneOid,
            seedBegin,
            seedEnd,
            uint32_t(seedEnd - seedBegin)});

        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        for(uint32_t position=seedBegin; position<seedEnd; position++) {
            const Shasta2AnchorId anchorId = backboneJourney[position];
            if(anchorOwner[uint64_t(anchorId)] == anchorUnclaimed) {
                anchorOwner[uint64_t(anchorId)] = windowId;
                ++task.claimedAnchorCount;
                ++backboneClaimedAnchors;
            }
        }

        ++epoch;
        touchedOrientedReads.clear();
        for(uint32_t position=seedBegin; position<seedEnd; position++) {
            const Shasta2AnchorId anchorId = backboneJourney[position];
            const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& ami: anchor) {
                const OrientedReadId oid = ami.orientedReadId;
                if(oid == backboneOid || ami.positionInJourney == invalid<uint32_t>) {
                    continue;
                }
                const uint32_t oidValue = uint32_t(oid.getValue());
                if(touchedEpoch[oidValue] != epoch) {
                    touchedEpoch[oidValue] = epoch;
                    touchedMin[oidValue] = ami.positionInJourney;
                    touchedMax[oidValue] = ami.positionInJourney;
                    touchedCount[oidValue] = 1;
                    touchedOrientedReads.push_back(oidValue);
                } else {
                    touchedMin[oidValue] = min(touchedMin[oidValue], ami.positionInJourney);
                    touchedMax[oidValue] = max(touchedMax[oidValue], ami.positionInJourney);
                    ++touchedCount[oidValue];
                }
            }
        }
        anchorWindowTouchedReads += touchedOrientedReads.size();

        for(const uint32_t oidValue: touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(oid.getValue() >= shasta2Journeys->size()) {
                continue;
            }
            const auto journey = (*shasta2Journeys)[oid];
            if(journey.empty()) {
                continue;
            }
            const uint32_t begin = touchedMin[oidValue];
            const uint32_t end = min<uint32_t>(touchedMax[oidValue] + 1, uint32_t(journey.size()));
            uint32_t position = begin;
            while(position < end) {
                while(position < end && anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                    ++position;
                }
                const uint32_t runBegin = position;
                while(position < end && anchorOwner[uint64_t(journey[position])] == anchorUnclaimed) {
                    anchorOwner[uint64_t(journey[position])] = windowId;
                    ++task.claimedAnchorCount;
                    ++nonBackboneClaimedAnchors;
                    ++position;
                }
                if(runBegin != position) {
                    task.readIntervals.push_back(AnchorWindowReadInterval{
                        oid,
                        runBegin,
                        position,
                        touchedCount[oidValue]});
                    ++anchorWindowSplitIntervals;
                }
            }
        }

        anchorWindowClaimedAnchors += task.claimedAnchorCount;
        anchorWindowReadIntervals += task.readIntervals.size();
        maxAnchorWindowReadIntervals = max<uint64_t>(
            maxAnchorWindowReadIntervals,
            task.readIntervals.size());
        maxAnchorWindowClaimedAnchors = max<uint64_t>(
            maxAnchorWindowClaimedAnchors,
            task.claimedAnchorCount);
        anchorWindows.push_back(std::move(task));
    };

    // Initialize the heap: all strand-0 reads, or a single diagnostic backbone only.
    if(restrictAnchorWindowPlannerToSingleBackboneReadId) {
        const OrientedReadId backboneOid(anchorWindowPlannerOnlyBackboneReadId, 0);
        if(backboneOid.getValue() >= shasta2Journeys->size()) {
            ++anchorWindowSkippedNoJourney;
        } else {
            const auto journey = (*shasta2Journeys)[backboneOid];
            if(journey.empty()) {
                ++anchorWindowSkippedNoJourney;
            } else {
                pushCandidate(backboneOid, journey, 0, uint32_t(journey.size()));
            }
        }
    } else {
        for(const ReadId readId: anchorReadsByLength) {
            const OrientedReadId backboneOid(readId, 0);
            if(backboneOid.getValue() >= shasta2Journeys->size()) {
                ++anchorWindowSkippedNoJourney;
                continue;
            }
            const auto journey = (*shasta2Journeys)[backboneOid];
            if(journey.empty()) {
                ++anchorWindowSkippedNoJourney;
                continue;
            }
            pushCandidate(backboneOid, journey, 0, uint32_t(journey.size()));
        }
    }

    while(!candidateHeap.empty()) {
        const AnchorWindowCandidate candidate = candidateHeap.top();
        candidateHeap.pop();
        ++candidateIntervalsPopped;

        const ReadId readId = candidate.backboneOrientedReadId.getReadId();
        if(candidate.generation != candidateGeneration[uint64_t(readId)]) {
            ++discardedOldGenerationCandidates;
            continue;
        }
        if(candidate.backboneOrientedReadId.getValue() >= shasta2Journeys->size()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }
        const auto journey = (*shasta2Journeys)[candidate.backboneOrientedReadId];
        if(journey.empty() || candidate.end > journey.size()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }

        bool isStillUnclaimed = true;
        for(uint32_t position=candidate.begin; position<candidate.end; position++) {
            if(anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                isStillUnclaimed = false;
                break;
            }
        }

        if(!isStillUnclaimed) {
            ++staleCandidateIntervals;
            ++candidateGeneration[uint64_t(readId)];
            pushCurrentUnclaimedIntervals(candidate.backboneOrientedReadId);
            continue;
        }

        ++anchorWindowBackboneIntervals;
        createAnchorWindow(candidate.backboneOrientedReadId, candidate.begin, candidate.end);
    }

    uint64_t unclaimedAnchorCount = 0;
    for(const uint32_t owner: anchorOwner) {
        if(owner == anchorUnclaimed) {
            ++unclaimedAnchorCount;
        }
    }

    array<uint64_t, histogramUpperBounds.size()> readIntervalsPerWindowHistogram = {};
    array<uint64_t, histogramUpperBounds.size()> claimedAnchorsPerWindowHistogram = {};
    array<uint64_t, histogramUpperBounds.size()> backboneAnchorSpanHistogram = {};
    for(const AnchorWindowTask& task: anchorWindows) {
        addToHistogram(readIntervalsPerWindowHistogram, task.readIntervals.size());
        addToHistogram(claimedAnchorsPerWindowHistogram, task.claimedAnchorCount);
        addToHistogram(backboneAnchorSpanHistogram, task.backboneEnd - task.backboneBegin);
    }

    cout << timestamp << "[TheseusReadWindowMSA] Top 10 read journey anchor usage:" << endl;
    for(uint64_t i=0; i<min<uint64_t>(10, anchorReadsByLength.size()); i++) {
        const ReadId readId = anchorReadsByLength[i];
        const OrientedReadId oid(readId, 0);
        uint64_t journeyAnchorCount = 0;
        uint64_t claimedJourneyAnchors = 0;
        uint64_t journeyBaseSpan = 0;
        uint64_t claimedRuns = 0;
        uint64_t unclaimedRuns = 0;
        if(oid.getValue() < shasta2Journeys->size()) {
            const auto journey = (*shasta2Journeys)[oid];
            journeyAnchorCount = journey.size();
            if(!journey.empty()) {
                journeyBaseSpan = intervalBaseSpan(oid, journey, 0, uint32_t(journey.size()));
            }
            bool previousClaimed = false;
            bool previousUnclaimed = false;
            for(const Shasta2AnchorId anchorId: journey) {
                const bool isClaimed = anchorOwner[uint64_t(anchorId)] != anchorUnclaimed;
                if(isClaimed) {
                    ++claimedJourneyAnchors;
                    if(!previousClaimed) {
                        ++claimedRuns;
                    }
                } else if(!previousUnclaimed) {
                    ++unclaimedRuns;
                }
                previousClaimed = isClaimed;
                previousUnclaimed = !isClaimed;
            }
        }
        cout << timestamp << "  rank=" << i
             << " readId=" << readId
             << " length=" << reads->getRead(readId).baseCount
             << " journeyAnchors=" << journeyAnchorCount
             << " journeyBaseSpan=" << journeyBaseSpan
             << " claimedJourneyAnchors=" << claimedJourneyAnchors
             << " unclaimedJourneyAnchors=" << (journeyAnchorCount - claimedJourneyAnchors)
             << " claimedFraction=" << (journeyAnchorCount == 0 ? 0. : double(claimedJourneyAnchors) / double(journeyAnchorCount))
             << " claimedRuns=" << claimedRuns
             << " unclaimedRuns=" << unclaimedRuns
             << endl;
    }

    const auto anchorWindowEnd = chrono::steady_clock::now();
    const double anchorWindowSeconds = chrono::duration<double>(anchorWindowEnd - anchorWindowBegin).count();
    const auto totalEndAnchorPrototype = chrono::steady_clock::now();
    const double totalAnchorPrototypeSeconds = chrono::duration<double>(totalEndAnchorPrototype - totalBegin).count();

    cout << timestamp << "[TheseusReadWindowMSA] Anchor-window prototype ends."
         << " reads=" << readCount
         << " anchors=" << anchorCount
         << " windows=" << anchorWindows.size()
         << " minBackboneWindowAnchors=" << minBackboneWindowAnchors
         << " claimedAnchors=" << anchorWindowClaimedAnchors
         << " backboneClaimedAnchors=" << backboneClaimedAnchors
         << " nonBackboneClaimedAnchors=" << nonBackboneClaimedAnchors
         << " unclaimedAnchors=" << unclaimedAnchorCount
         << " backboneIntervals=" << anchorWindowBackboneIntervals
         << " candidateIntervalsPushed=" << candidateIntervalsPushed
         << " candidateIntervalsPopped=" << candidateIntervalsPopped
         << " staleCandidateIntervals=" << staleCandidateIntervals
         << " discardedOldGenerationCandidates=" << discardedOldGenerationCandidates
         << " readIntervals=" << anchorWindowReadIntervals
         << " splitReadIntervals=" << anchorWindowSplitIntervals
         << " touchedReads=" << anchorWindowTouchedReads
         << " maxReadIntervalsPerWindow=" << maxAnchorWindowReadIntervals
         << " avgReadIntervalsPerWindow=" << (anchorWindows.empty() ? 0. : double(anchorWindowReadIntervals) / double(anchorWindows.size()))
         << " maxClaimedAnchorsPerWindow=" << maxAnchorWindowClaimedAnchors
         << " avgClaimedAnchorsPerWindow=" << (anchorWindows.empty() ? 0. : double(anchorWindowClaimedAnchors) / double(anchorWindows.size()))
         << " readIntervalsPerWindowHistogram=" << histogramToString(readIntervalsPerWindowHistogram)
         << " claimedAnchorsPerWindowHistogram=" << histogramToString(claimedAnchorsPerWindowHistogram)
         << " backboneAnchorSpanHistogram=" << histogramToString(backboneAnchorSpanHistogram)
         << " skippedNoJourney=" << anchorWindowSkippedNoJourney
         << " skippedShortRuns=" << anchorWindowSkippedShortRuns
         << " anchorWindowSeconds=" << fixed << setprecision(6) << anchorWindowSeconds
         << " totalSeconds=" << totalAnchorPrototypeSeconds
         << defaultfloat << endl;

    if(runAnchorWindowBackbonePairTheseusMsa) {
        BackbonePairMsaCounters pairMsaCounters;
        const auto pairMsaBegin = chrono::steady_clock::now();
        runShasta2BackbonePairTheseusMsas(
            *this,
            anchorWindows,
            *shasta2Anchors,
            *shasta2Journeys,
            threadCount,
            orientedReadCount,
            pairMsaCounters);
        const auto pairMsaEnd = chrono::steady_clock::now();
        const double pairWallSeconds =
            chrono::duration<double>(pairMsaEnd - pairMsaBegin).count();
        cout << timestamp << "[TheseusReadWindowMSA] Backbone-pair Theseus MSA ends."
             << " pairJobsScheduled=" << pairMsaCounters.pairJobsScheduled
             << " msasRun=" << pairMsaCounters.msasRun
             << " skippedCoverage=" << pairMsaCounters.skippedCoverage
             << " skippedOrdinal=" << pairMsaCounters.skippedOrdinal
             << " skippedEmptyFocalWindow=" << pairMsaCounters.skippedEmptyFocalWindow
             << " skippedSingleSequence=" << pairMsaCounters.skippedSingleSequence
             << " skippedReverseOrderBothAnchorReads="
             << pairMsaCounters.skippedReverseOrderBothAnchorReads
             << " leftOnlySequences=" << pairMsaCounters.leftOnlySequences
             << " rightOnlySequences=" << pairMsaCounters.rightOnlySequences
             << " failedLeftOnlyAnchorChecks=" << pairMsaCounters.failedLeftOnlyAnchorChecks
             << " failedRightOnlyAnchorChecks=" << pairMsaCounters.failedRightOnlyAnchorChecks
             << " msaCpuSeconds=" << fixed << setprecision(6) << pairMsaCounters.msaCpuSeconds
             << " wallSeconds=" << pairWallSeconds
             << defaultfloat
             << endl;
    }
    return;
    }

    vector<ReadId> readsByLength;
    readsByLength.reserve(readCount);
    for(uint64_t readId=0; readId<readCount; readId++) {
        readsByLength.push_back(ReadId(readId));
    }
    sort(readsByLength.begin(), readsByLength.end(),
        [&](ReadId a, ReadId b) {
            const uint64_t lengthA = reads->getRead(a).baseCount;
            const uint64_t lengthB = reads->getRead(b).baseCount;
            if(lengthA != lengthB) {
                return lengthA > lengthB;
            }
            return a < b;
        });

    const uint32_t unclaimed = numeric_limits<uint32_t>::max();
    vector<uint32_t> readOwner(readCount, unclaimed);
    vector<ReadWindowTask> windows;
    uint64_t crossWindowEdgeCount = 0;
    uint64_t claimedReadCount = 0;
    uint64_t scannedReadGraphEdges = 0;
    uint64_t skippedCrossStrandEdges = 0;
    uint64_t skippedInconsistentEdges = 0;
    uint64_t skippedSelfEdges = 0;
    uint64_t borrowedReadCount = 0;
    uint64_t rejectedBackboneCandidates = 0;
    uint64_t rejectedClaimedNeighborEdges = 0;

    const auto planBegin = chrono::steady_clock::now();
    for(const ReadId seedReadId: readsByLength) {
        if(readOwner[uint64_t(seedReadId)] != unclaimed) {
            continue;
        }

        vector<CandidateEvidence> candidateEvidence;
        bool touchesClaimedRead = false;
        const OrientedReadId seedOid(seedReadId, 0);
        for(const uint32_t edgeId: readGraph.connectivity[seedOid.getValue()]) {
            ++scannedReadGraphEdges;
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            if(edge.crossesStrands) {
                ++skippedCrossStrandEdges;
                continue;
            }
            if(edge.hasInconsistentAlignment) {
                ++skippedInconsistentEdges;
                continue;
            }
            const OrientedReadId other = edge.getOther(seedOid);
            const uint64_t otherReadId = uint64_t(other.getReadId());
            if(otherReadId == uint64_t(seedReadId)) {
                ++skippedSelfEdges;
                continue;
            }

            if(readOwner[otherReadId] != unclaimed) {
                touchesClaimedRead = true;
                ++rejectedClaimedNeighborEdges;
            }
            candidateEvidence.push_back(CandidateEvidence{
                other,
                uint32_t(edge.alignmentId)});
        }

        if(touchesClaimedRead) {
            ++rejectedBackboneCandidates;
            continue;
        }

        const uint32_t windowId = uint32_t(windows.size());
        ReadWindowTask task;
        task.windowId = windowId;
        task.backboneReadId = seedReadId;
        task.orientedReads.push_back(OrientedReadId(seedReadId, 0));
        task.claimedReads.push_back(seedReadId);
        readOwner[uint64_t(seedReadId)] = windowId;
        ++claimedReadCount;

        for(const CandidateEvidence& evidence: candidateEvidence) {
            task.orientedReads.push_back(evidence.orientedReadId);
            task.alignmentIds.push_back(evidence.alignmentId);

            const ReadId otherReadId = evidence.orientedReadId.getReadId();
            if(readOwner[uint64_t(otherReadId)] == unclaimed) {
                readOwner[uint64_t(otherReadId)] = windowId;
                ++claimedReadCount;
                task.claimedReads.push_back(otherReadId);
            }
        }

        windows.push_back(std::move(task));
    }
    const auto planEnd = chrono::steady_clock::now();
    const double planSeconds = chrono::duration<double>(planEnd - planBegin).count();

    const auto anchorRescueBegin = chrono::steady_clock::now();
    vector<uint32_t> sharedAnchorCount(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint64_t anchorRescueRows = 0;
    uint64_t anchorRescueWindows = 0;
    uint64_t anchorRescueSharedAnchorHits = 0;
    uint64_t anchorRescueSkippedNoJourney = 0;
    uint64_t anchorRescueSkippedClaimedReads = 0;
    uint64_t anchorRescueSkippedLowSharedAnchors = 0;

    for(ReadWindowTask& task: windows) {
        const OrientedReadId backboneOid(task.backboneReadId, 0);
        if(backboneOid.getValue() >= shasta2Journeys->size()) {
            ++anchorRescueSkippedNoJourney;
            continue;
        }

        const auto journey = (*shasta2Journeys)[backboneOid];
        if(journey.empty()) {
            ++anchorRescueSkippedNoJourney;
            continue;
        }

        touchedOrientedReads.clear();
        for(const Shasta2AnchorId anchorId: journey) {
            const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& ami: anchor) {
                const OrientedReadId oid = ami.orientedReadId;
                if(oid == backboneOid) {
                    continue;
                }
                const ReadId readId = oid.getReadId();
                if(readOwner[uint64_t(readId)] != unclaimed) {
                    continue;
                }
                const uint64_t oidValue = uint64_t(oid.getValue());
                if(sharedAnchorCount[oidValue] == 0) {
                    touchedOrientedReads.push_back(uint32_t(oidValue));
                }
                ++sharedAnchorCount[oidValue];
                ++anchorRescueSharedAnchorHits;
            }
        }

        uint64_t rescuedInWindow = 0;
        for(const uint32_t oidValue: touchedOrientedReads) {
            const uint32_t shared = sharedAnchorCount[oidValue];
            sharedAnchorCount[oidValue] = 0;
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(readOwner[uint64_t(oid.getReadId())] != unclaimed) {
                ++anchorRescueSkippedClaimedReads;
                continue;
            }
            if(shared < minSharedAnchorsForRescue) {
                ++anchorRescueSkippedLowSharedAnchors;
                continue;
            }

            task.orientedReads.push_back(oid);
            task.alignmentIds.push_back(invalidAlignmentId);
            ++anchorRescueRows;
            ++rescuedInWindow;
        }
        if(rescuedInWindow != 0) {
            ++anchorRescueWindows;
        }
    }
    const auto anchorRescueEnd = chrono::steady_clock::now();
    const double anchorRescueSeconds = chrono::duration<double>(anchorRescueEnd - anchorRescueBegin).count();

    uint64_t singletonWindowCount = 0;
    uint64_t maxClaimedReadCount = 0;
    uint64_t totalEvidenceReadCount = 0;
    uint64_t maxEvidenceReadCount = 0;
    uint64_t ownerMismatchCount = 0;
    vector<bool> isBackboneRead(readCount, false);
    for(const ReadWindowTask& task: windows) {
        maxClaimedReadCount = max<uint64_t>(maxClaimedReadCount, task.claimedReads.size());
        maxEvidenceReadCount = max<uint64_t>(maxEvidenceReadCount, task.orientedReads.size());
        totalEvidenceReadCount += task.orientedReads.size();
        if(task.claimedReads.size() == 1) {
            ++singletonWindowCount;
        }
        isBackboneRead[uint64_t(task.backboneReadId)] = true;
        for(const ReadId readId: task.claimedReads) {
            if(readOwner[uint64_t(readId)] != task.windowId) {
                ++ownerMismatchCount;
            }
        }
    }

    uint64_t backboneConflictEdgeCount = 0;
    for(const ReadWindowTask& task: windows) {
        const OrientedReadId seedOid(task.backboneReadId, 0);
        for(const uint32_t edgeId: readGraph.connectivity[seedOid.getValue()]) {
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                continue;
            }
            const ReadId otherReadId = edge.getOther(seedOid).getReadId();
            if(isBackboneRead[uint64_t(otherReadId)] &&
                readOwner[uint64_t(otherReadId)] != task.windowId) {
                ++backboneConflictEdgeCount;
            }
        }
    }
    backboneConflictEdgeCount /= 2;

    vector<vector<EvidenceOccurrence> > orientedReadOccurrences(2 * readCount);
    for(const ReadWindowTask& task: windows) {
        for(uint64_t row=0; row<task.orientedReads.size(); row++) {
            const OrientedReadId oid = task.orientedReads[row];
            const uint32_t alignmentId = (row == 0) ?
                invalidAlignmentId :
                task.alignmentIds[row - 1];
            orientedReadOccurrences[uint64_t(oid.getValue())].push_back(EvidenceOccurrence{
                task.windowId,
                uint32_t(row),
                alignmentId});
        }
    }

    uint64_t indexedOrientedReadCount = 0;
    uint64_t sharedOrientedReadCount = 0;
    uint64_t maxOrientedReadOccurrenceCount = 0;
    for(const vector<EvidenceOccurrence>& occurrences: orientedReadOccurrences) {
        if(occurrences.empty()) {
            continue;
        }
        ++indexedOrientedReadCount;
        maxOrientedReadOccurrenceCount = max<uint64_t>(
            maxOrientedReadOccurrenceCount,
            occurrences.size());
        if(occurrences.size() > 1) {
            ++sharedOrientedReadCount;
        }
    }

    const ReadId debugReadId = ReadId(2109);
    const uint32_t debugOwner = readOwner[uint64_t(debugReadId)];
    cout << timestamp << "[TheseusReadWindowMSA] Debug physical read " << debugReadId << "\n"
         << timestamp << "  length=" << reads->getRead(debugReadId).baseCount << "\n"
         << timestamp << "  ownerWindow=" << debugOwner;
    if(debugOwner != unclaimed) {
        const ReadWindowTask& ownerTask = windows[debugOwner];
        cout << "\n"
             << timestamp << "  ownerBackbone=" << OrientedReadId(ownerTask.backboneReadId, 0) << "\n"
             << timestamp << "  ownerBackboneLength=" << reads->getRead(ownerTask.backboneReadId).baseCount << "\n"
             << timestamp << "  ownerClaimedReads=" << ownerTask.claimedReads.size() << "\n"
             << timestamp << "  ownerEvidenceRows=" << ownerTask.orientedReads.size();
    }
    cout << endl;

    vector<bool> printedDebugWindow(windows.size(), false);
    uint64_t debugOccurrenceCount = 0;
    for(Strand strand=0; strand<2; strand++) {
        const OrientedReadId debugOrientedReadId(debugReadId, strand);
        const vector<EvidenceOccurrence>& debugOccurrences =
            orientedReadOccurrences[uint64_t(debugOrientedReadId.getValue())];
        cout << timestamp << "[TheseusReadWindowMSA] Occurrences of "
             << debugOrientedReadId
             << " count=" << debugOccurrences.size();
        for(const EvidenceOccurrence& occurrence: debugOccurrences) {
            ++debugOccurrenceCount;
            cout << "\n"
                 << timestamp << "  windowId=" << occurrence.windowId
                 << " row=" << occurrence.row;
            if(occurrence.alignmentId == invalidAlignmentId) {
            cout << " noAlignmentId";
            } else {
                cout << " alignmentId=" << occurrence.alignmentId;
            }
        }
        cout << endl;
    }

    for(Strand strand=0; strand<2; strand++) {
        const OrientedReadId debugOrientedReadId(debugReadId, strand);
        const vector<EvidenceOccurrence>& debugOccurrences =
            orientedReadOccurrences[uint64_t(debugOrientedReadId.getValue())];
        for(const EvidenceOccurrence& occurrence: debugOccurrences) {
            if(printedDebugWindow[occurrence.windowId]) {
                continue;
            }
            printedDebugWindow[occurrence.windowId] = true;

            const ReadWindowTask& task = windows[occurrence.windowId];
            cout << timestamp << "[TheseusReadWindowMSA] Window containing physical read "
                 << debugReadId << "\n"
                 << timestamp << "  windowId=" << task.windowId << "\n"
                 << timestamp << "  backbone=" << OrientedReadId(task.backboneReadId, 0) << "\n"
                 << timestamp << "  backboneLength=" << reads->getRead(task.backboneReadId).baseCount << "\n"
                 << timestamp << "  claimedReads=" << task.claimedReads.size() << "\n"
                 << timestamp << "  evidenceRows=" << task.orientedReads.size() << "\n"
                 << timestamp << "  debugReadOwner=" << debugOwner << "\n"
                 << timestamp << "  isOwnerWindow=" << (task.windowId == debugOwner)
                 << endl;

            printWrappedItems(cout, "[TheseusReadWindowMSA] Claimed reads:", task.claimedReads, 24);

            cout << timestamp << "[TheseusReadWindowMSA] Evidence rows:";
            for(uint64_t row=0; row<task.orientedReads.size(); row++) {
                if((row % 8) == 0) {
                    cout << "\n" << timestamp << "  ";
                }
                cout << " " << row << ":" << task.orientedReads[row];
                if(row > 0) {
                    if(task.alignmentIds[row - 1] == invalidAlignmentId) {
                        cout << "(anchorRescue)";
                    } else {
                        cout << "(alignmentId=" << task.alignmentIds[row - 1] << ")";
                    }
                } else {
                    cout << "(backbone)";
                }
            }
            cout << endl;
        }
    }
    if(debugOccurrenceCount == 0) {
        cout << timestamp << "[TheseusReadWindowMSA] No window contains physical read "
             << debugReadId
             << " as an evidence row."
             << endl;
    }

    vector<ThreadCounters> threadCounters(threadCount);
    double msaWallSeconds = 0.;
#if 0
    atomic<uint64_t> nextWindow(0);
    const auto msaBegin = chrono::steady_clock::now();
    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        threads.emplace_back([&, threadId]() {
            ThreadCounters& counters = threadCounters[threadId];
            theseus::Penalties penalties(0, 2, 3, 1);
            theseus::Heuristics heuristics(false, false);
            while(true) {
                const uint64_t windowIndex = nextWindow.fetch_add(1);
                if(windowIndex >= windows.size()) {
                    break;
                }
                const ReadWindowTask& task = windows[windowIndex];
                if(task.orientedReads.size() < 2) {
                    ++counters.skippedSmallWindows;
                    continue;
                }

                vector<string> sequences;
                sequences.reserve(task.orientedReads.size());
                uint64_t baseCount = 0;
                for(const OrientedReadId oid: task.orientedReads) {
                    sequences.push_back(extractWholeOrientedReadSequence(getReads(), oid));
                    baseCount += sequences.back().size();
                }
                if(sequences.empty() || sequences.front().empty()) {
                    ++counters.skippedSmallWindows;
                    continue;
                }

                const auto begin = chrono::steady_clock::now();
                theseus::TheseusMSA aligner(
                    penalties,
                    heuristics,
                    sequences.front(),
                    1,
                    false);
                for(size_t i=1; i<sequences.size(); i++) {
                    if(!sequences[i].empty()) {
                        aligner.align(sequences[i], 1, false, true);
                    }
                }
                ostringstream discard;
                aligner.print_as_msa(discard);
                const auto end = chrono::steady_clock::now();

                ++counters.windows;
                counters.rows += sequences.size();
                counters.bases += baseCount;
                counters.msaSeconds += chrono::duration<double>(end - begin).count();
            }
        });
    }
    for(thread& t: threads) {
        t.join();
    }
    const auto msaEnd = chrono::steady_clock::now();
    msaWallSeconds = chrono::duration<double>(msaEnd - msaBegin).count();
#endif

    ThreadCounters totalCounters;
    for(const ThreadCounters& counters: threadCounters) {
        totalCounters.windows += counters.windows;
        totalCounters.skippedSmallWindows += counters.skippedSmallWindows;
        totalCounters.rows += counters.rows;
        totalCounters.bases += counters.bases;
        totalCounters.msaSeconds += counters.msaSeconds;
    }

    const auto totalEnd = chrono::steady_clock::now();
    const double totalSeconds = chrono::duration<double>(totalEnd - totalBegin).count();
    cout << timestamp << "[TheseusReadWindowMSA] Prototype ends."
         << " reads=" << readCount
         << " alignments=" << alignmentData.size()
         << " readGraphEdges=" << readGraph.edges.size()
         << " scannedReadGraphEdges=" << scannedReadGraphEdges
         << " skippedCrossStrandEdges=" << skippedCrossStrandEdges
         << " skippedInconsistentEdges=" << skippedInconsistentEdges
         << " skippedSelfEdges=" << skippedSelfEdges
         << " windows=" << windows.size()
         << " claimedReads=" << claimedReadCount
         << " unclaimedReads=" << (readCount - claimedReadCount)
         << " rejectedBackboneCandidates=" << rejectedBackboneCandidates
         << " rejectedClaimedNeighborEdges=" << rejectedClaimedNeighborEdges
         << " anchorRescueMinSharedAnchors=" << minSharedAnchorsForRescue
         << " anchorRescueRows=" << anchorRescueRows
         << " anchorRescueWindows=" << anchorRescueWindows
         << " anchorRescueSharedAnchorHits=" << anchorRescueSharedAnchorHits
         << " anchorRescueSkippedNoJourney=" << anchorRescueSkippedNoJourney
         << " anchorRescueSkippedClaimedReads=" << anchorRescueSkippedClaimedReads
         << " anchorRescueSkippedLowSharedAnchors=" << anchorRescueSkippedLowSharedAnchors
         << " singletonWindows=" << singletonWindowCount
         << " maxClaimedReadsPerWindow=" << maxClaimedReadCount
         << " avgClaimedReadsPerWindow=" << (windows.empty() ? 0. : double(claimedReadCount) / double(windows.size()))
         << " evidenceReads=" << totalEvidenceReadCount
         << " borrowedEvidenceReads=" << borrowedReadCount
         << " maxEvidenceReadsPerWindow=" << maxEvidenceReadCount
         << " avgEvidenceReadsPerWindow=" << (windows.empty() ? 0. : double(totalEvidenceReadCount) / double(windows.size()))
         << " indexedOrientedReads=" << indexedOrientedReadCount
         << " sharedOrientedReads=" << sharedOrientedReadCount
         << " maxOrientedReadOccurrences=" << maxOrientedReadOccurrenceCount
         << " crossWindowEdges=" << crossWindowEdgeCount
         << " backboneConflictEdges=" << backboneConflictEdgeCount
         << " ownerMismatches=" << ownerMismatchCount
         << " runMsa=" << runTheseusMsa
         << " processedWindows=" << totalCounters.windows
         << " skippedSmallWindows=" << totalCounters.skippedSmallWindows
         << " rows=" << totalCounters.rows
         << " bases=" << totalCounters.bases
         << " planSeconds=" << fixed << setprecision(6) << planSeconds
         << " anchorRescueSeconds=" << anchorRescueSeconds
         << " msaThreadSeconds=" << totalCounters.msaSeconds
         << " msaWallSeconds=" << msaWallSeconds
         << " totalSeconds=" << totalSeconds
         << " avgSecondsPerWindow=" << (windows.empty() ? 0. : totalSeconds / double(windows.size()))
         << " avgMsaThreadSecondsPerWindow=" << (totalCounters.windows ? totalCounters.msaSeconds / double(totalCounters.windows) : 0.)
         << " threadCount=" << threadCount
         << defaultfloat << endl;
}
