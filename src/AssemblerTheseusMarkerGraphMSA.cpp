// AssemblerTheseusMarkerGraphMSA.cpp
//
// Diagnostic prototype: for oriented read 0-0, run Theseus/Pericles MSAs over
// consecutive marker-graph vertex anchors and print variation sites detected in
// the MSA. This intentionally does not feed EC/parity or AlignedEvidenceStore.

#include "Assembler.hpp"
#include "Base.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "Reads.hpp"
#include "findMarkerId.hpp"
#include "mode3-Anchor.hpp"
#include "timestamp.hpp"

#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace dinara;
using namespace dinara::mode3;
using namespace std;

namespace {

constexpr uint64_t oneSidedOffsetRatioNumerator = 11;
constexpr uint64_t oneSidedOffsetRatioDenominator = 10;
constexpr uint64_t minAnchorCoverageForMsa = 6;
constexpr uint64_t minSnpRefSupport = 3;
constexpr uint64_t minSnpAltSupport = 3;
constexpr uint64_t minReportedAltLength = 16;
constexpr uint64_t minFilteredHomopolymerRunLength = 3;
constexpr bool writeSiteMsaFiles = false;
constexpr bool logSites = false;
constexpr bool logPairSummaries = false;
constexpr bool logRecruitmentAudit = false;
constexpr bool logDebugMsa = false;
const string siteMsaOutputDirectory = "TheseusMGMSA-sites";
constexpr AnchorId debugMsaLeftAnchor = 64014;
constexpr AnchorId debugMsaRightAnchor = 64016;
constexpr uint32_t auditRecruitmentTargetPosition = 43191;

struct SequenceInfo {
    OrientedReadId oid;
    string sequence;
    uint32_t begin = 0;
    uint32_t end = 0;
    bool hasBothAnchors = false;
    bool isEndsFree = false;
    char anchorSide = 'B'; // B=both, L=left only, R=right only.
};

struct Segment {
    string sequence;
    uint32_t begin = 0;
    uint32_t end = 0;
};

struct AnchorPairKey {
    AnchorId left = 0;
    AnchorId right = 0;

    bool operator==(const AnchorPairKey& that) const
    {
        return left == that.left && right == that.right;
    }
};

struct AnchorPairKeyHash {
    size_t operator()(const AnchorPairKey& key) const
    {
        return (std::hash<uint64_t>{}(uint64_t(key.left)) << 1) ^
            std::hash<uint64_t>{}(uint64_t(key.right));
    }
};

string displayAllele(const string& allele)
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

bool isCanonicalBase(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

vector<string> parseMsaFasta(const string& text)
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

string orientedReadList(const vector<uint64_t>& indexes, const vector<SequenceInfo>& sequenceInfos)
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

string orientedReadValueList(const vector<uint64_t>& values)
{
    string s;
    for(size_t i=0; i<values.size(); i++) {
        if(i) {
            s += ",";
        }
        s += OrientedReadId::fromValue(ReadId(values[i])).getString();
    }
    return s;
}

void addReadIndex(vector<uint64_t>& indexes, uint64_t index)
{
    if(find(indexes.begin(), indexes.end(), index) == indexes.end()) {
        indexes.push_back(index);
    }
}

string alleleType(const string& ref, const string& alt)
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

bool shouldReportAllele(const string& type, const string& alt)
{
    return type == "SNP" || alt.size() >= minReportedAltLength;
}

bool isHomopolymerAt(const string& seq, uint32_t pos)
{
    if(pos >= seq.size()) {
        return false;
    }

    const size_t p = size_t(pos);
    const char base = seq[p];
    if(!isCanonicalBase(base)) {
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
    return end - begin >= minFilteredHomopolymerRunLength;
}

bool snpTouchesHomopolymer(const string& focalTargetSequence, uint32_t pos, char altBase)
{
    if(isHomopolymerAt(focalTargetSequence, pos)) {
        return true;
    }
    if(pos >= focalTargetSequence.size() || !isCanonicalBase(altBase)) {
        return false;
    }
    string altTargetSequence = focalTargetSequence;
    altTargetSequence[pos] = altBase;
    return isHomopolymerAt(altTargetSequence, pos);
}

string siteMsaFileName(
    uint64_t pairIndex,
    AnchorId left,
    AnchorId right,
    uint32_t targetPos)
{
    ostringstream name;
    name << siteMsaOutputDirectory
         << "/site-pair" << pairIndex
         << "-anchors" << left << "-" << right
         << "-pos" << targetPos
         << ".msa.fa";
    return name.str();
}

void printVariationSitesFromMsa(
    uint64_t pairIndex,
    AnchorId left,
    AnchorId right,
    const vector<SequenceInfo>& sequenceInfos,
    const vector<string>& alignedSequences,
    const string& msaText,
    uint32_t focalBegin,
    uint32_t focalEnd,
    uint32_t reportBegin,
    uint32_t reportEnd,
    double msaSeconds)
{
    if(alignedSequences.size() != sequenceInfos.size() || alignedSequences.empty()) {
        cout << timestamp << "[TheseusMGMSA] pair=" << pairIndex
             << " skipped: MSA sequence count mismatch." << endl;
        return;
    }

    const size_t columnCount = alignedSequences.front().size();
    for(const string& s: alignedSequences) {
        if(s.size() != columnCount) {
            cout << timestamp << "[TheseusMGMSA] pair=" << pairIndex
                 << " skipped: MSA rows have inconsistent lengths." << endl;
            return;
        }
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
        focalTargetSequence.push_back(isCanonicalBase(base) ? base : 'N');
    }

    vector<uint8_t> isDirtyTargetOffset(columnByTargetOffset.size(), 0);
    for(size_t offset=0; offset<columnByTargetOffset.size(); offset++) {
        const size_t c = columnByTargetOffset[offset];
        const char refBase = alignedSequences[0][c];
        if(!isCanonicalBase(refBase)) {
            continue;
        }
        for(size_t r=1; r<alignedSequences.size(); r++) {
            if(!hasBase[r] || c < firstNonGap[r] || c > lastNonGap[r]) {
                continue;
            }
            const char base = alignedSequences[r][c];
            if(base == '-') {
                isDirtyTargetOffset[offset] = 1;
            } else if(isCanonicalBase(base) && base != refBase) {
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
                        if(isCanonicalBase(alignedSequences[r][j])) {
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
    uint64_t eventCount = 0;
    map<string, uint64_t> typeCounts;
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
            if(isCanonicalBase(base)) {
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
                if(isCanonicalBase(base)) {
                    allele.push_back(base);
                }
            }
            if(allele == ref) {
                refReads.push_back(r);
            } else {
                const AlleleKey key{alleleType(ref, allele), ref, allele};
                addReadIndex(altReadsByAllele[key], r);
            }
        }

        uint64_t maxAltSupport = 0;
        map<AlleleKey, vector<uint64_t>> reportableAltReadsByAllele;
        for(const auto& [key, altReads]: altReadsByAllele) {
            const auto& [type, ref, alt] = key;
            if(!shouldReportAllele(type, alt)) {
                continue;
            }
            if(type == "SNP" &&
               snpTouchesHomopolymer(focalTargetSequence, uint32_t(beginOffset), alt[0])) {
                continue;
            }
            reportableAltReadsByAllele[key] = altReads;
            maxAltSupport = max<uint64_t>(maxAltSupport, altReads.size());
        }
        if(refReads.size() < minSnpRefSupport || maxAltSupport < minSnpAltSupport) {
            beginOffset = endOffset + 1;
            continue;
        }

        if(logSites) {
            cout << "[TheseusMGMSA] SITE"
                 << " pair=" << pairIndex
                 << " anchors=" << left << "->" << right
                 << " focal=" << sequenceInfos.front().oid
                 << " pos=" << targetPos
                 << " ref=" << displayAllele(ref)
                 << " refN=" << refReads.size()
                 << " refReads=" << orientedReadList(refReads, sequenceInfos);
            for(const auto& [key, altReads]: reportableAltReadsByAllele) {
                const auto& [type, ref, alt] = key;
                cout << " alt=" << type << ":" << displayAllele(ref) << ">"
                     << displayAllele(alt)
                     << ":N=" << altReads.size()
                     << ":reads=" << orientedReadList(altReads, sequenceInfos);
            }
            cout << endl;
        }
        for(const auto& [key, altReads]: reportableAltReadsByAllele) {
            (void) altReads;
            const auto& [type, ref, alt] = key;
            ++eventCount;
            ++typeCounts[type];
        }

        if(writeSiteMsaFiles) {
            filesystem::create_directories(siteMsaOutputDirectory);
            const string fileName = siteMsaFileName(pairIndex, left, right, targetPos);
            const string tmpFileName = fileName + ".tmp";
            {
                ofstream file(tmpFileName, ios::binary);
                file << msaText;
                file.close();
            }
            error_code ec;
            filesystem::rename(tmpFileName, fileName, ec);
            if(!ec) {
                cout << "[TheseusMGMSA] SITE_MSA"
                     << " pair=" << pairIndex
                     << " anchors=" << left << "->" << right
                     << " pos=" << targetPos
                     << " file=" << fileName
                     << endl;
            } else {
                cout << timestamp << "[TheseusMGMSA] failed to write site MSA file "
                     << fileName
                     << " error=" << ec.message()
                     << endl;
            }
        }
        beginOffset = endOffset + 1;
    }

    if(logPairSummaries) {
        cout << timestamp << "[TheseusMGMSA] pair=" << pairIndex
             << " anchors=" << left << "->" << right
             << " focalRange=" << focalBegin << "-" << focalEnd
             << " sequences=" << sequenceInfos.size()
             << " events=" << eventCount
             << " seconds=" << fixed << setprecision(6) << msaSeconds;
        for(const auto& [type, count]: typeCounts) {
            cout << " " << type << "=" << count;
        }
        cout << defaultfloat << endl;
    }
}

Segment extractSegmentFromOrdinals(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    Segment segment;
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
        return Segment{};
    }

    segment.sequence.reserve(segment.end - segment.begin);
    for(uint32_t pos=segment.begin; pos<segment.end; pos++) {
        segment.sequence.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return segment;
}

Segment extractSegmentFromBases(
    const Reads& reads,
    OrientedReadId oid,
    uint32_t begin,
    uint32_t end)
{
    Segment segment;
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

bool segmentContainsMarkerAtLeft(
    const Segment& segment,
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

bool segmentContainsMarkerAtRight(
    const Segment& segment,
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

struct SharedAnchorOnRead {
    uint32_t targetJourneyIndex = 0;
    AnchorId anchorId = 0;
    uint32_t targetOrdinal = 0;
    uint32_t readOrdinal = 0;
};

pair<uint32_t, uint32_t> bestIncreasingRun(span<const SharedAnchorOnRead> sharedAnchors)
{
    uint32_t bestBegin = invalid<uint32_t>;
    uint32_t bestEnd = invalid<uint32_t>;
    uint32_t runBegin = 0;

    for(uint32_t i=1; i<=sharedAnchors.size(); i++) {
        const bool continues = (i < sharedAnchors.size()) &&
            (sharedAnchors[i - 1].readOrdinal < sharedAnchors[i].readOrdinal);
        if(continues) {
            continue;
        }

        if(i - runBegin >= 2) {
            if(bestBegin == invalid<uint32_t> ||
               i - runBegin > bestEnd - bestBegin ||
               (i - runBegin == bestEnd - bestBegin &&
                sharedAnchors[i - 1].targetOrdinal - sharedAnchors[runBegin].targetOrdinal >
                sharedAnchors[bestEnd - 1].targetOrdinal - sharedAnchors[bestBegin].targetOrdinal)) {
                bestBegin = runBegin;
                bestEnd = i;
            }
        }
        runBegin = i;
    }

    return {bestBegin, bestEnd};
}

} // namespace



void Assembler::computeTheseusMarkerGraphMSAPrototype(
    uint64_t maxAnchorPairs,
    uint64_t maxReadsPerPair,
    uint64_t threadCount)
{
    cout << timestamp << "[TheseusMGMSA] Prototype begins for all oriented reads." << endl;
    const auto prototypeBegin = chrono::steady_clock::now();

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    DINARA_ASSERT(assemblerInfo.isOpen);
    DINARA_ASSERT((assemblerInfo->k % 2) == 0);
    DINARA_ASSERT(reads->readCount() > 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    auto anchors = make_shared<Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        2,
        numeric_limits<uint64_t>::max(),
        threadCount,
        true);
    anchors->computeJourneys(threadCount);

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    unordered_set<AnchorPairKey, AnchorPairKeyHash> processedPairs;

    uint64_t processedPairCount = 0;
    uint64_t processedReadCount = 0;
    uint64_t skippedNoJourneyStorage = 0;
    uint64_t skippedShortJourney = 0;
    uint64_t skippedShortFilteredJourney = 0;
    uint64_t totalEventBearingPairs = 0;
    uint64_t skippedSingleSequencePairs = 0;
    uint64_t skippedReverseOrderBothAnchorReads = 0;
    uint64_t leftOnlySequenceCount = 0;
    uint64_t rightOnlySequenceCount = 0;
    uint64_t failedLeftOnlyAnchorChecks = 0;
    uint64_t failedRightOnlyAnchorChecks = 0;

    for(uint64_t readId=0; readId<readCount && processedPairCount<maxAnchorPairs; readId++) {
        const OrientedReadId focalOid(ReadId(readId), 0);
        if(focalOid.getValue() >= anchors->journeys.size()) {
            ++skippedNoJourneyStorage;
            continue;
        }

        const auto focalJourney = anchors->journeys[focalOid.getValue()];
        if(focalJourney.size() < 2) {
            ++skippedShortJourney;
            continue;
        }

        vector<AnchorId> filteredFocalJourney;
        filteredFocalJourney.reserve(focalJourney.size());
        for(uint32_t i=0; i<focalJourney.size(); i++) {
            const AnchorId anchorId = focalJourney[i];
            const Anchor anchor = (*anchors)[anchorId];
            if(anchor.size() >= minAnchorCoverageForMsa) {
                filteredFocalJourney.push_back(anchorId);
            }
        }
        if(filteredFocalJourney.size() < 2) {
            ++skippedShortFilteredJourney;
            continue;
        }
        ++processedReadCount;

        for(uint32_t focalJourneyPos=0;
            focalJourneyPos + 1 < filteredFocalJourney.size() && processedPairCount < maxAnchorPairs;
            focalJourneyPos++) {

        const AnchorId left = filteredFocalJourney[focalJourneyPos];
        const AnchorId right = filteredFocalJourney[focalJourneyPos + 1];
        const AnchorPairKey key{left, right};
        if(!processedPairs.insert(key).second) {
            continue;
        }

        const uint32_t focalLeftOrdinal = anchors->getFirstOrdinal(left, focalOid);
        const uint32_t focalRightOrdinal = anchors->getFirstOrdinal(right, focalOid);
        const Anchor leftAnchor = (*anchors)[left];
        const Anchor rightAnchor = (*anchors)[right];

        Segment focalWindowSegment = extractSegmentFromOrdinals(
            getReads(), *markers, assemblerInfo->k, focalOid, focalLeftOrdinal, focalRightOrdinal);
        if(focalWindowSegment.sequence.empty()) {
            continue;
        }

        unordered_map<uint64_t, uint32_t> leftOrdinals;
        unordered_map<uint64_t, uint32_t> rightOrdinals;
        for(const AnchorMarkerInterval& interval: leftAnchor) {
            leftOrdinals.try_emplace(interval.orientedReadId.getValue(), interval.ordinal0);
        }
        for(const AnchorMarkerInterval& interval: rightAnchor) {
            rightOrdinals.try_emplace(interval.orientedReadId.getValue(), interval.ordinal0);
        }

        vector<uint64_t> anchorUnionValues;
        anchorUnionValues.reserve(leftOrdinals.size() + rightOrdinals.size());
        for(const auto& [oidValue, ordinal]: leftOrdinals) {
            (void) ordinal;
            anchorUnionValues.push_back(oidValue);
        }
        for(const auto& [oidValue, ordinal]: rightOrdinals) {
            (void) ordinal;
            if(!leftOrdinals.contains(oidValue)) {
                anchorUnionValues.push_back(oidValue);
            }
        }
        sort(anchorUnionValues.begin(), anchorUnionValues.end());

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
            const auto readMarkers = (*markers)[oid.getValue()];
            const uint32_t rightOrdinal = itRight->second;
            if(leftOrdinal >= readMarkers.size() || rightOrdinal >= readMarkers.size()) {
                continue;
            }
            if(leftOrdinal >= rightOrdinal) {
                continue;
            }
            const uint32_t kHalf = uint32_t(assemblerInfo->k / 2);
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
        Segment focalSegment = extractSegmentFromBases(
            getReads(),
            focalOid,
            (focalWindowSegment.begin > alignmentPadding) ?
                (focalWindowSegment.begin - alignmentPadding) : 0,
            focalWindowSegment.end + alignmentPadding);
        if(focalSegment.sequence.empty()) {
            continue;
        }

        vector<SequenceInfo> sequenceInfos;
        sequenceInfos.reserve(size_t(min<uint64_t>(maxReadsPerPair, 1024)));
        sequenceInfos.push_back(SequenceInfo{
            focalOid,
            focalSegment.sequence,
            focalSegment.begin,
            focalSegment.end,
            true,
            false,
            'B'});

        uint64_t auditSkippedReverseOrder = 0;
        uint64_t auditSkippedEmpty = 0;
        uint64_t auditSkippedAnchorCheck = 0;
        uint64_t auditSkippedOutOfRange = 0;
        for(const uint64_t oidValue64: candidateValues) {
            if(sequenceInfos.size() >= maxReadsPerPair) {
                break;
            }
            if(oidValue64 >= orientedReadCount) {
                ++auditSkippedOutOfRange;
                continue;
            }

            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue64));
            const bool hasLeft = leftOrdinals.contains(oidValue64);
            const bool hasRight = rightOrdinals.contains(oidValue64);
            Segment segment;
            bool hasBothAnchors = false;
            bool isEndsFree = true;
            char anchorSide = 'B';

            if(hasLeft && hasRight) {
                if(leftOrdinals[oidValue64] >= rightOrdinals[oidValue64]) {
                    ++skippedReverseOrderBothAnchorReads;
                    ++auditSkippedReverseOrder;
                    continue;
                }
                const auto readMarkers = (*markers)[oid.getValue()];
                const uint32_t leftOrdinal = leftOrdinals[oidValue64];
                const uint32_t rightOrdinal = rightOrdinals[oidValue64];
                if(leftOrdinal < readMarkers.size() && rightOrdinal < readMarkers.size()) {
                    const uint32_t kHalf = uint32_t(assemblerInfo->k / 2);
                    const uint32_t windowBegin = readMarkers[leftOrdinal].position + kHalf;
                    const uint32_t windowEnd = readMarkers[rightOrdinal].position + kHalf;
                    if(windowEnd > windowBegin) {
                        segment = extractSegmentFromBases(
                            getReads(),
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
                const auto readMarkers = (*markers)[oid.getValue()];
                const uint32_t leftOrdinal = leftOrdinals[oidValue64];
                if(leftOrdinal < readMarkers.size()) {
                    const uint32_t begin =
                        readMarkers[leftOrdinal].position + uint32_t(assemblerInfo->k / 2);
                    segment = extractSegmentFromBases(
                        getReads(), oid, begin, begin + approximateSpan);
                }
                anchorSide = 'L';
            }

            if(segment.sequence.empty() && hasRight) {
                const uint32_t rightOrdinal = rightOrdinals[oidValue64];
                const auto readMarkers = (*markers)[oid.getValue()];
                if(rightOrdinal < readMarkers.size()) {
                    const uint32_t end =
                        readMarkers[rightOrdinal].position + uint32_t(assemblerInfo->k / 2);
                    const uint32_t begin = (end > approximateSpan) ? (end - approximateSpan) : 0;
                    segment = extractSegmentFromBases(getReads(), oid, begin, end);
                }
                anchorSide = 'R';
            }

            if(segment.sequence.empty()) {
                ++auditSkippedEmpty;
                continue;
            }

            if(anchorSide == 'L') {
                if(!segmentContainsMarkerAtLeft(
                    segment, *markers, assemblerInfo->k, oid, leftOrdinals[oidValue64])) {
                    ++failedLeftOnlyAnchorChecks;
                    ++auditSkippedAnchorCheck;
                    continue;
                }
                ++leftOnlySequenceCount;
            } else if(anchorSide == 'R') {
                if(!segmentContainsMarkerAtRight(
                    segment, *markers, assemblerInfo->k, oid, rightOrdinals[oidValue64])) {
                    ++failedRightOnlyAnchorChecks;
                    ++auditSkippedAnchorCheck;
                    continue;
                }
                ++rightOnlySequenceCount;
            }

            sequenceInfos.push_back(SequenceInfo{
                oid,
                segment.sequence,
                segment.begin,
                segment.end,
                hasBothAnchors,
                isEndsFree,
                anchorSide});
        }

        if(logRecruitmentAudit &&
           focalWindowSegment.begin <= auditRecruitmentTargetPosition &&
           auditRecruitmentTargetPosition < focalWindowSegment.end) {
            unordered_set<uint64_t> addedValues;
            for(const SequenceInfo& info: sequenceInfos) {
                addedValues.insert(info.oid.getValue());
            }
            vector<uint64_t> missingValues;
            for(const uint64_t oidValue: anchorUnionValues) {
                if(!addedValues.contains(oidValue)) {
                    missingValues.push_back(oidValue);
                }
            }
            vector<uint64_t> addedSorted(addedValues.begin(), addedValues.end());
            sort(addedSorted.begin(), addedSorted.end());

            cout << "[TheseusMGMSA] RECRUITMENT_AUDIT"
                 << " pair=" << processedPairCount
                 << " anchors=" << left << "->" << right
                 << " targetPos=" << auditRecruitmentTargetPosition
                 << " focalRange=" << focalWindowSegment.begin << "-" << focalWindowSegment.end
                 << " leftReads=" << leftOrdinals.size()
                 << " rightReads=" << rightOrdinals.size()
                 << " unionReads=" << anchorUnionValues.size()
                 << " sequenceInfos=" << sequenceInfos.size()
                 << " maxReadsPerPair=" << maxReadsPerPair
                 << " skippedReverseOrder=" << auditSkippedReverseOrder
                 << " skippedEmpty=" << auditSkippedEmpty
                 << " skippedAnchorCheck=" << auditSkippedAnchorCheck
                 << " skippedOutOfRange=" << auditSkippedOutOfRange
                 << " missingFromMsa=" << missingValues.size()
                 << endl;
            cout << "[TheseusMGMSA] RECRUITMENT_AUDIT_UNION"
                 << " anchors=" << left << "->" << right
                 << " targetPos=" << auditRecruitmentTargetPosition
                 << " reads=" << orientedReadValueList(anchorUnionValues)
                 << endl;
            cout << "[TheseusMGMSA] RECRUITMENT_AUDIT_ADDED"
                 << " anchors=" << left << "->" << right
                 << " targetPos=" << auditRecruitmentTargetPosition
                 << " reads=" << orientedReadValueList(addedSorted)
                 << endl;
            if(!missingValues.empty()) {
                cout << "[TheseusMGMSA] RECRUITMENT_AUDIT_MISSING"
                     << " anchors=" << left << "->" << right
                     << " targetPos=" << auditRecruitmentTargetPosition
                     << " reads=" << orientedReadValueList(missingValues)
                     << endl;
            }
        }

        if(sequenceInfos.size() < 2) {
            ++skippedSingleSequencePairs;
            continue;
        }

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
        const string msaText = msaOut.str();
        const auto msaEnd = chrono::steady_clock::now();
        const double msaSeconds = chrono::duration<double>(msaEnd - msaBegin).count();

        const vector<string> alignedSequences = parseMsaFasta(msaText);
        if(logDebugMsa && left == debugMsaLeftAnchor && right == debugMsaRightAnchor) {
            cout << "[TheseusMGMSA] DEBUG_MSA_BEGIN"
                 << " pair=" << processedPairCount
                 << " anchors=" << left << "->" << right
                 << " focalRange=" << focalWindowSegment.begin << "-" << focalWindowSegment.end
                 << " sequences=" << sequenceInfos.size()
                 << " alignedSequences=" << alignedSequences.size()
                 << endl;
            for(size_t i=0; i<sequenceInfos.size(); i++) {
                const SequenceInfo& info = sequenceInfos[i];
                cout << ">seq" << i
                     << " read=" << info.oid
                     << " range=" << info.begin << "-" << info.end
                     << " anchors=" << (info.hasBothAnchors ? "2end" : "1end")
                     << " side=" << info.anchorSide
                     << " endsFree=" << (info.isEndsFree ? 1 : 0)
                     << endl;
                if(i < alignedSequences.size()) {
                    cout << alignedSequences[i] << endl;
                } else {
                    cout << "<missing aligned sequence>" << endl;
                }
            }
            cout << "[TheseusMGMSA] DEBUG_MSA_END"
                 << " anchors=" << left << "->" << right
                 << endl;
        }
        printVariationSitesFromMsa(
            processedPairCount,
            left,
            right,
            sequenceInfos,
            alignedSequences,
            msaText,
            focalSegment.begin,
            focalSegment.end,
            focalWindowSegment.begin,
            focalWindowSegment.end,
            msaSeconds);

        ++totalEventBearingPairs;
        ++processedPairCount;
    }
    }

    const auto prototypeEnd = chrono::steady_clock::now();
    const double prototypeSeconds = chrono::duration<double>(prototypeEnd - prototypeBegin).count();
    cout << timestamp << "[TheseusMGMSA] Prototype ends."
         << " reads=" << readCount
         << " processedReads=" << processedReadCount
         << " skippedNoJourneyStorage=" << skippedNoJourneyStorage
         << " skippedShortJourney=" << skippedShortJourney
         << " skippedShortFilteredJourney=" << skippedShortFilteredJourney
         << " processedPairs=" << processedPairCount
         << " skippedSingleSequencePairs=" << skippedSingleSequencePairs
         << " skippedReverseOrderBothAnchorReads=" << skippedReverseOrderBothAnchorReads
         << " leftOnlySequences=" << leftOnlySequenceCount
         << " rightOnlySequences=" << rightOnlySequenceCount
         << " failedLeftOnlyAnchorChecks=" << failedLeftOnlyAnchorChecks
         << " failedRightOnlyAnchorChecks=" << failedRightOnlyAnchorChecks
         << " reportedPairs=" << totalEventBearingPairs
         << " seconds=" << fixed << setprecision(6) << prototypeSeconds
         << " avgSecondsPerRead=" << (readCount ? prototypeSeconds / double(readCount) : 0.)
         << " avgSecondsPerProcessedRead=" << (processedReadCount ? prototypeSeconds / double(processedReadCount) : 0.)
         << defaultfloat << endl;
}



void Assembler::computeTheseusTargetBackboneMSAPrototype(
    uint64_t maxReads,
    uint64_t threadCount)
{
    cout << timestamp << "[TheseusBackboneMSA] Prototype begins for oriented read 0-0." << endl;

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    DINARA_ASSERT(assemblerInfo.isOpen);
    DINARA_ASSERT((assemblerInfo->k % 2) == 0);
    DINARA_ASSERT(reads->readCount() > 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    auto anchors = make_shared<Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        2,
        numeric_limits<uint64_t>::max(),
        threadCount,
        true);
    anchors->computeJourneys(threadCount);

    const OrientedReadId focalOid(ReadId(0), 0);
    if(focalOid.getValue() >= anchors->journeys.size()) {
        cout << timestamp << "[TheseusBackboneMSA] oriented read 0-0 has no journey storage." << endl;
        return;
    }

    const auto focalJourney = anchors->journeys[focalOid.getValue()];
    vector<AnchorId> filteredFocalJourney;
    vector<uint32_t> filteredFocalOrdinals;
    filteredFocalJourney.reserve(focalJourney.size());
    filteredFocalOrdinals.reserve(focalJourney.size());
    for(uint32_t i=0; i<focalJourney.size(); i++) {
        const AnchorId anchorId = focalJourney[i];
        const Anchor anchor = (*anchors)[anchorId];
        if(anchor.size() >= minAnchorCoverageForMsa) {
            filteredFocalJourney.push_back(anchorId);
            filteredFocalOrdinals.push_back(anchors->getFirstOrdinal(anchorId, focalOid));
        }
    }
    cout << timestamp << "[TheseusBackboneMSA] read 0-0 journey length=" << focalJourney.size()
         << " filteredJourneyLength=" << filteredFocalJourney.size()
         << " minAnchorCoverage=" << minAnchorCoverageForMsa << endl;
    if(filteredFocalJourney.size() < 2) {
        cout << timestamp << "[TheseusBackboneMSA] read 0-0 filtered journey too short." << endl;
        return;
    }

    unordered_map<uint64_t, vector<SharedAnchorOnRead>> sharedAnchorsByRead;
    for(uint32_t targetIndex=0; targetIndex<filteredFocalJourney.size(); targetIndex++) {
        const AnchorId anchorId = filteredFocalJourney[targetIndex];
        const Anchor anchor = (*anchors)[anchorId];
        for(const AnchorMarkerInterval& interval: anchor) {
            const uint64_t oidValue = interval.orientedReadId.getValue();
            if(oidValue == focalOid.getValue()) {
                continue;
            }
            sharedAnchorsByRead[oidValue].push_back(SharedAnchorOnRead{
                targetIndex,
                anchorId,
                filteredFocalOrdinals[targetIndex],
                interval.ordinal0});
        }
    }

    const uint32_t targetReadLength =
        uint32_t(getReads().getRead(focalOid.getReadId()).baseCount);
    Segment focalSegment = extractSegmentFromBases(getReads(), focalOid, 0, targetReadLength);
    if(focalSegment.sequence.empty()) {
        cout << timestamp << "[TheseusBackboneMSA] focal backbone sequence is empty." << endl;
        return;
    }

    vector<SequenceInfo> sequenceInfos;
    sequenceInfos.reserve(size_t(min<uint64_t>(maxReads, sharedAnchorsByRead.size() + 1)));
    sequenceInfos.push_back(SequenceInfo{
        focalOid,
        focalSegment.sequence,
        focalSegment.begin,
        focalSegment.end,
        true,
        false,
        'B'});

    vector<uint64_t> readValues;
    readValues.reserve(sharedAnchorsByRead.size());
    for(const auto& [oidValue, sharedAnchors]: sharedAnchorsByRead) {
        (void) sharedAnchors;
        readValues.push_back(oidValue);
    }
    sort(readValues.begin(), readValues.end());

    uint64_t skippedTooFewSharedAnchors = 0;
    uint64_t skippedNoIncreasingRun = 0;
    uint64_t skippedEmptySegment = 0;
    uint64_t skippedOutOfRange = 0;
    uint64_t addedReadSegments = 0;
    uint64_t sharedAnchorCount = 0;
    uint64_t chosenSharedAnchorCount = 0;
    uint64_t chosenTargetSpanSum = 0;
    uint64_t chosenReadSpanSum = 0;
    AnchorId firstBackboneAnchor = filteredFocalJourney.front();
    AnchorId lastBackboneAnchor = filteredFocalJourney.back();

    const uint64_t orientedReadCount = 2 * reads->readCount();
    for(const uint64_t oidValue: readValues) {
        if(sequenceInfos.size() >= maxReads) {
            break;
        }
        if(oidValue >= orientedReadCount) {
            ++skippedOutOfRange;
            continue;
        }

        vector<SharedAnchorOnRead>& sharedAnchors = sharedAnchorsByRead[oidValue];
        sharedAnchorCount += sharedAnchors.size();
        if(sharedAnchors.size() < 2) {
            ++skippedTooFewSharedAnchors;
            continue;
        }

        sort(sharedAnchors.begin(), sharedAnchors.end(),
            [](const SharedAnchorOnRead& a, const SharedAnchorOnRead& b) {
                if(a.targetJourneyIndex != b.targetJourneyIndex) {
                    return a.targetJourneyIndex < b.targetJourneyIndex;
                }
                return a.readOrdinal < b.readOrdinal;
            });
        const auto [runBegin, runEnd] = bestIncreasingRun(span<const SharedAnchorOnRead>(
            sharedAnchors.data(), sharedAnchors.size()));
        if(runBegin == invalid<uint32_t>) {
            ++skippedNoIncreasingRun;
            continue;
        }

        const SharedAnchorOnRead& first = sharedAnchors[runBegin];
        const SharedAnchorOnRead& last = sharedAnchors[runEnd - 1];
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        Segment segment = extractSegmentFromOrdinals(
            getReads(), *markers, assemblerInfo->k, oid, first.readOrdinal, last.readOrdinal);
        if(segment.sequence.empty()) {
            ++skippedEmptySegment;
            continue;
        }

        chosenSharedAnchorCount += runEnd - runBegin;
        chosenTargetSpanSum += (last.targetOrdinal > first.targetOrdinal) ?
            (last.targetOrdinal - first.targetOrdinal) : 0;
        chosenReadSpanSum += (last.readOrdinal > first.readOrdinal) ?
            (last.readOrdinal - first.readOrdinal) : 0;
        firstBackboneAnchor = min(firstBackboneAnchor, first.anchorId);
        lastBackboneAnchor = max(lastBackboneAnchor, last.anchorId);

        sequenceInfos.push_back(SequenceInfo{
            oid,
            segment.sequence,
            segment.begin,
            segment.end,
            true,
            true,
            'B'});
        ++addedReadSegments;
    }

    if(sequenceInfos.size() < 2) {
        cout << timestamp << "[TheseusBackboneMSA] skipped: only focal sequence available."
             << " candidateReads=" << sharedAnchorsByRead.size()
             << " skippedTooFewSharedAnchors=" << skippedTooFewSharedAnchors
             << " skippedNoIncreasingRun=" << skippedNoIncreasingRun
             << " skippedEmptySegment=" << skippedEmptySegment
             << endl;
        return;
    }

    cout << timestamp << "[TheseusBackboneMSA] sequences=" << sequenceInfos.size()
         << " addedReadSegments=" << addedReadSegments
         << " candidateReads=" << sharedAnchorsByRead.size()
         << " sharedAnchors=" << sharedAnchorCount
         << " chosenSharedAnchors=" << chosenSharedAnchorCount
         << " chosenTargetOrdinalSpanSum=" << chosenTargetSpanSum
         << " chosenReadOrdinalSpanSum=" << chosenReadSpanSum
         << " skippedTooFewSharedAnchors=" << skippedTooFewSharedAnchors
         << " skippedNoIncreasingRun=" << skippedNoIncreasingRun
         << " skippedEmptySegment=" << skippedEmptySegment
         << " skippedOutOfRange=" << skippedOutOfRange
         << endl;

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
        aligner.align(sequenceInfos[i].sequence, 1, false, true);
    }

    ostringstream msaOut;
    aligner.print_as_msa(msaOut);
    const string msaText = msaOut.str();
    const auto msaEnd = chrono::steady_clock::now();
    const double msaSeconds = chrono::duration<double>(msaEnd - msaBegin).count();
    const vector<string> alignedSequences = parseMsaFasta(msaText);

    printVariationSitesFromMsa(
        0,
        firstBackboneAnchor,
        lastBackboneAnchor,
        sequenceInfos,
        alignedSequences,
        msaText,
        focalSegment.begin,
        focalSegment.end,
        focalSegment.begin,
        focalSegment.end,
        msaSeconds);

    cout << timestamp << "[TheseusBackboneMSA] Prototype ends."
         << " sequences=" << sequenceInfos.size()
         << " msaSeconds=" << fixed << setprecision(6) << msaSeconds
         << defaultfloat << endl;
}
