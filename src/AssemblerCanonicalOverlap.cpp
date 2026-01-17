// AssemblerCanonicalOverlap.cpp
// Functions for building and using the canonical per-Read overlap index.

#include "Assembler.hpp"
#include "CanonicalOverlap.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
#include "AssemblerPhasing.hpp"
#include <thread>

using namespace dinara;

// Build the CanonicalOverlap index from alignmentData.
// This converts the per-OrientedRead AlignmentData to per-Read CanonicalOverlap.
void Assembler::buildCanonicalOverlapIndex()
{
    cout << timestamp << "Building canonical overlap index." << endl;
    
    const uint32_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();
    
    // Initialize the index
    overlapIndex.resize(readCount);
    
    // Convert each AlignmentData to CanonicalOverlap
    for (uint64_t i = 0; i < alignmentCount; i++) {
        const auto& ad = alignmentData[i];
        
        CanonicalOverlap ov;
        
        // Read IDs (AlignmentData guarantees readIds[0] < readIds[1])
        ov.qn = ad.readIds[0];
        ov.tn = ad.readIds[1];
        
        // Coordinates
        ov.qs = ad.qs;
        ov.qe = ad.qe;
        ov.ts = ad.ts;
        ov.te = ad.te;
        
        // Strand (isSameStrand in AlignmentData)
        ov.rev = ad.isSameStrand ? 0 : 1;
        
        // Marker info from AlignmentInfo
        ov.markerCount = (uint16_t)std::min(ad.info.markerCount, (uint32_t)65535);
        ov.firstOrdinal0 = (uint16_t)std::min(ad.info.data[0].firstOrdinal, (uint32_t)65535);
        ov.lastOrdinal0 = (uint16_t)std::min(ad.info.data[0].lastOrdinal, (uint32_t)65535);
        ov.firstOrdinal1 = (uint16_t)std::min(ad.info.data[1].firstOrdinal, (uint32_t)65535);
        
        // Flags: Initialize to not deleted, unknown phasing
        ov.del = 0;
        ov.is_match = 0;  // Unknown
        ov.strong = 0;
        ov.lg_indel = ad.hasLargeIndel ? 1 : 0;
        
        // Chain count / alignment length
        ov.cc = ad.info.markerCount;
        
        // CIGAR index (same as alignment index for now)
        ov.cigarIdx = (uint32_t)i;
        
        // Add to index (ensures canonical ordering qn < tn)
        overlapIndex.addOverlap(ov);
    }
    
    cout << timestamp << "Built canonical overlap index with " 
         << overlapIndex.totalOverlaps << " overlaps for " 
         << readCount << " reads." << endl;
}


// Perform phasing using CanonicalOverlap storage (sets is_match/strong flags)
// HIFIASM PARITY: Uses score-based classification (not kept/not-kept)
void Assembler::performPhasingCanonical(uint64_t threadCount)
{
    cout << timestamp << "Performing canonical phasing (is_match/strong flags)." << endl;
    
    checkPhasingCigarsAreOpen();
    
    const uint32_t readCount = reads->readCount();
    
    // Phasing config - use struct defaults for ONT mode
    PhasingConfig config;
    config.hom_cov = assemblerInfo->kmerDistributionInfo.coveragePeak * 2;
    config.is_ont = true;
    // Note: Use default values from PhasingConfig for other parameters:
    // s_hap_cov=3, infor_cov=1, st_rate=0.05, st_max=2
    
    // Statistics
    uint64_t totalCis = 0, totalTrans = 0, totalStrong = 0, totalWeak = 0;
    
    for (uint32_t readId = 0; readId < readCount; readId++) {
        auto& overlaps = overlapIndex.getOverlapsAsQuery(readId);
        
        if (overlaps.empty()) continue;
        
        // Build PhasingOverlap list for non-deleted overlaps
        std::vector<PhasingOverlap> phasingOverlaps;
        std::vector<uint32_t> overlapIndices;  // Map from phasingOverlaps index to overlaps index
        phasingOverlaps.reserve(overlaps.size());
        overlapIndices.reserve(overlaps.size());
        
        for (uint32_t idx = 0; idx < overlaps.size(); idx++) {
            const auto& ov = overlaps[idx];
            if (ov.isDeleted()) continue;
            
            PhasingOverlap po;
            po.alnIdx = ov.cigarIdx;
            po.targetReadId = readId;
            po.queryReadId = ov.tn;
            po.queryStrand = ov.rev;
            po.targetStart = ov.qs;
            po.targetEnd = ov.qe;
            po.queryStart = ov.ts;
            po.queryEnd = ov.te;
            
            auto cigarSpan = phasingCigars[ov.cigarIdx];
            po.cigar.assign(cigarSpan.begin(), cigarSpan.end());
            
            phasingOverlaps.push_back(std::move(po));
            overlapIndices.push_back(idx);
        }
        
        if (phasingOverlaps.empty()) continue;
        
        // Get phasing scores for ALL overlaps (Hifiasm parity)
        auto phasingResults = AssemblerPhasing::getOverlapPhasingScores(
            *this, readId, phasingOverlaps, config
        );
        
        // Apply scores to CanonicalOverlap is_match and strong flags
        for (size_t i = 0; i < phasingResults.size(); i++) {
            uint32_t ovIdx = overlapIndices[i];
            auto& ov = overlaps[ovIdx];
            
            const auto& result = phasingResults[i];
            
            // Hifiasm parity:
            // score > 0 → TRANS (is_match = 2)
            // score <= 0 → CIS (is_match = 1)
            if (result.score > 0) {
                ov.is_match = 2;  // TRANS
                totalTrans++;
            } else {
                ov.is_match = 1;  // CIS
                totalCis++;
            }
            
            // strong = hasInformativeSite (Hifiasm parity)
            if (result.hasInformativeSite) {
                ov.strong = 1;
                totalStrong++;
            } else {
                ov.strong = 0;
                totalWeak++;
            }
        }
    }
    
    // DEBUG: Output CIS overlaps for read 0
    cout << timestamp << "DEBUG: Read 0 CIS overlaps:" << endl;
    auto& read0Overlaps = overlapIndex.getOverlapsAsQuery(0);
    for (const auto& ov : read0Overlaps) {
        if (ov.isCis() && !ov.isDeleted()) {
            // ov.rev: 0 = same strand, 1 = reverse
            // ov.cc = marker count (alignment length), NOT phasing score
            cout << "  Read 0-0 CIS with Read " << ov.tn << "-" << (ov.rev ? "1" : "0") 
                 << " (strong=" << (int)ov.strong << ", markers=" << ov.cc << ")" << endl;
        }
    }
    // Also check as target
    overlapIndex.forEachOverlapAsTarget(0, [](const CanonicalOverlap& ov) {
        if (ov.isCis() && !ov.isDeleted()) {
            cout << "  Read 0-0 CIS with Read " << ov.qn << "-" << (ov.rev ? "1" : "0")
                 << " (strong=" << (int)ov.strong << ")" << endl;
        }
    });
    
    cout << timestamp << "Canonical phasing complete: " 
         << totalCis << " CIS, " << totalTrans << " TRANS, "
         << totalStrong << " strong, " << totalWeak << " weak." << endl;
}

// Thread function for parallel phasing (placeholder for now)
void Assembler::performPhasingCanonicalThreadFunction(uint64_t threadId)
{
    // TODO: Implement parallel version if needed
    // For now, performPhasingCanonical is single-threaded
}


// Create read graph using canonical OverlapIndex storage
// Uses is_match for phasing (1=CIS, 2=TRANS) and del for filtering
void Assembler::createReadGraph7()
{
    cout << timestamp << "createReadGraph7 begins (using OverlapIndex)." << endl;
    
    const uint32_t readCount = reads->readCount();
    
    // Count statistics
    uint64_t totalOverlaps = overlapIndex.totalOverlaps;
    uint64_t transCount = 0;
    uint64_t deletedCount = 0;
    uint64_t palindromicCount = 0;
    uint64_t chimericCount = 0;
    uint64_t keptCount = 0;
    
    // Collect overlaps to include in read graph
    std::vector<bool> keepOverlap(totalOverlaps, false);
    
    // Traverse all overlaps
    uint64_t overlapIdx = 0;
    for (uint32_t qn = 0; qn < readCount; qn++) {
        for (auto& ov : overlapIndex.sources[qn]) {
            // 1. Check phasing: TRANS overlaps are excluded
            if (ov.isTrans()) {
                transCount++;
                overlapIdx++;
                continue;
            }
            
            // 2. Check deletion (from filtering)
            if (ov.isDeleted()) {
                deletedCount++;
                overlapIdx++;
                continue;
            }
            
            // 3. Check palindromic reads
            if (reads->getFlags(ov.qn).isPalindromic || 
                reads->getFlags(ov.tn).isPalindromic) {
                palindromicCount++;
                overlapIdx++;
                continue;
            }
            
            // 4. Check chimeric reads
            if (isChimericRead.size() > 0) {
                if (isChimericRead[ov.qn] || isChimericRead[ov.tn]) {
                    chimericCount++;
                    overlapIdx++;
                    continue;
                }
            }
            
            // Keep this overlap
            keepOverlap[overlapIdx] = true;
            keptCount++;
            overlapIdx++;
        }
    }
    
    cout << timestamp << "TRANS phasing removed " << transCount << " overlaps." << endl;
    cout << timestamp << "Filtering removed " << deletedCount << " overlaps." << endl;
    cout << timestamp << "Palindromic filter removed " << palindromicCount << " overlaps." << endl;
    cout << timestamp << "Chimeric filter removed " << chimericCount << " overlaps." << endl;
    cout << timestamp << "Kept " << keptCount << " / " << totalOverlaps << " overlaps." << endl;
    
    // Create read graph from kept overlaps
    // Convert back to alignmentData indices for compatibility
    std::vector<bool> keepAlignment(alignmentData.size(), false);
    
    overlapIdx = 0;
    for (uint32_t qn = 0; qn < readCount; qn++) {
        for (const auto& ov : overlapIndex.sources[qn]) {
            if (keepOverlap[overlapIdx]) {
                // Use cigarIdx to map back to alignmentData
                if (ov.cigarIdx < alignmentData.size()) {
                    keepAlignment[ov.cigarIdx] = true;
                    alignmentData[ov.cigarIdx].info.isInReadGraph = 1;
                }
            }
            overlapIdx++;
        }
    }
    
    createReadGraphUsingSelectedAlignments(keepAlignment);
    
    cout << timestamp << "createReadGraph7 completed." << endl;
}


// Create read graph directly from canonical OverlapIndex without mapping to alignmentData.
// This is Option A: Use overlapIndex directly for graph connectivity.
// Each canonical overlap (qn < tn) generates 2 edges in the OrientedRead graph.
void Assembler::createReadGraphFromCanonicalOverlaps()
{
    cout << timestamp << "createReadGraphFromCanonicalOverlaps begins." << endl;
    
    const uint32_t readCount = reads->readCount();
    
    // Count kept overlaps
    uint64_t keptCount = 0;
    uint64_t transCount = 0;
    uint64_t deletedCount = 0;
    uint64_t palindromicCount = 0;
    uint64_t chimericCount = 0;
    
    // First pass: count edges
    for (uint32_t qn = 0; qn < readCount; qn++) {
        for (const auto& ov : overlapIndex.sources[qn]) {
            if (ov.isTrans()) { transCount++; continue; }
            if (ov.isDeleted()) { deletedCount++; continue; }
            if (reads->getFlags(ov.qn).isPalindromic || 
                reads->getFlags(ov.tn).isPalindromic) { palindromicCount++; continue; }
            if (isChimericRead.size() > 0 &&
                (isChimericRead[ov.qn] || isChimericRead[ov.tn])) { chimericCount++; continue; }
            keptCount++;
        }
    }
    
    cout << timestamp << "Overlap filtering:" << endl;
    cout << "  TRANS removed: " << transCount << endl;
    cout << "  Deleted: " << deletedCount << endl;
    cout << "  Palindromic: " << palindromicCount << endl;
    cout << "  Chimeric: " << chimericCount << endl;
    cout << "  Kept: " << keptCount << endl;
    
    // Create read graph edges (2 edges per overlap: forward + RC)
    readGraph.edges.createNew(largeDataName("ReadGraphEdges"), largeDataPageSize);
    readGraph.edges.reserve(keptCount * 2);
    
    // Second pass: create edges
    for (uint32_t qn = 0; qn < readCount; qn++) {
        for (const auto& ov : overlapIndex.sources[qn]) {
            // Apply same filters
            if (ov.isTrans()) continue;
            if (ov.isDeleted()) continue;
            if (reads->getFlags(ov.qn).isPalindromic || 
                reads->getFlags(ov.tn).isPalindromic) continue;
            if (isChimericRead.size() > 0 &&
                (isChimericRead[ov.qn] || isChimericRead[ov.tn])) continue;
            
            // Create edge: (qn-strand0) to (tn-strand based on rev)
            ReadGraphEdge edge;
            edge.alignmentId = ov.cigarIdx & 0x3fff'ffff'ffff'ffff;  // Use cigarIdx as reference
            edge.crossesStrands = 0;
            edge.hasInconsistentAlignment = 0;
            edge.orientedReadIds[0] = OrientedReadId(ov.qn, 0);
            edge.orientedReadIds[1] = OrientedReadId(ov.tn, ov.rev ? 1 : 0);
            
            // Ensure canonical ordering for edge (orientedReadIds[0] < orientedReadIds[1])
            if (edge.orientedReadIds[0] > edge.orientedReadIds[1]) {
                std::swap(edge.orientedReadIds[0], edge.orientedReadIds[1]);
            }
            readGraph.edges.push_back(edge);
            
            // Also create the reverse complemented edge
            edge.orientedReadIds[0].flipStrand();
            edge.orientedReadIds[1].flipStrand();
            if (edge.orientedReadIds[0] > edge.orientedReadIds[1]) {
                std::swap(edge.orientedReadIds[0], edge.orientedReadIds[1]);
            }
            readGraph.edges.push_back(edge);
        }
    }
    
    readGraph.edges.unreserve();
    
    // Create read graph connectivity
    readGraph.connectivity.createNew(largeDataName("ReadGraphConnectivity"), largeDataPageSize);
    readGraph.connectivity.beginPass1(2 * readCount);
    for (const ReadGraphEdge& edge : readGraph.edges) {
        readGraph.connectivity.incrementCount(edge.orientedReadIds[0].getValue());
        readGraph.connectivity.incrementCount(edge.orientedReadIds[1].getValue());
    }
    readGraph.connectivity.beginPass2();
    for (size_t i = 0; i < readGraph.edges.size(); i++) {
        const ReadGraphEdge& edge = readGraph.edges[i];
        readGraph.connectivity.store(edge.orientedReadIds[0].getValue(), uint32_t(i));
        readGraph.connectivity.store(edge.orientedReadIds[1].getValue(), uint32_t(i));
    }
    readGraph.connectivity.endPass2();
    
    // Count isolated reads
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;
    for (ReadId readId = 0; readId < readCount; readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        if (readGraph.connectivity.size(orientedReadId.getValue()) == 0) {
            isolatedReadCount++;
            isolatedReadBaseCount += reads->getReadRawSequenceLength(readId);
        }
    }
    assemblerInfo->isolatedReadCount = isolatedReadCount;
    assemblerInfo->isolatedReadBaseCount = isolatedReadBaseCount;
    
    cout << timestamp << "Created read graph with " << readGraph.edges.size() << " edges." << endl;
    cout << timestamp << "Isolated reads: " << isolatedReadCount << endl;
    cout << timestamp << "createReadGraphFromCanonicalOverlaps completed." << endl;
}


// Create read graph from all filtered alignments (after filterBestHitAlignments).
// Simply keeps all alignments where isDeleted() is false.
void Assembler::createReadGraphFromFilteredAlignments()
{
    cout << timestamp << "createReadGraphFromFilteredAlignments begins." << endl;
    
    const uint32_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();
    
    // Count kept alignments
    uint64_t keptCount = 0;
    uint64_t deletedCount = 0;
    uint64_t palindromicCount = 0;
    
    // First pass: count edges
    for (uint64_t i = 0; i < alignmentCount; i++) {
        const auto& ad = alignmentData[i];
        
        if (ad.isDeleted()) { deletedCount++; continue; }
        if (reads->getFlags(ad.readIds[0]).isPalindromic || 
            reads->getFlags(ad.readIds[1]).isPalindromic) { palindromicCount++; continue; }
        
        keptCount++;
    }
    
    cout << timestamp << "Alignment filtering:" << endl;
    cout << "  Deleted by filterBestHitAlignments: " << deletedCount << endl;
    cout << "  Palindromic: " << palindromicCount << endl;
    cout << "  Kept: " << keptCount << " / " << alignmentCount << endl;
    
    // Create read graph edges (2 edges per alignment: forward + RC)
    readGraph.edges.createNew(largeDataName("ReadGraphEdges"), largeDataPageSize);
    readGraph.edges.reserve(keptCount * 2);
    
    // Second pass: create edges
    for (uint64_t i = 0; i < alignmentCount; i++) {
        auto& ad = alignmentData[i];
        
        // Skip deleted/palindromic
        if (ad.isDeleted()) continue;
        if (reads->getFlags(ad.readIds[0]).isPalindromic || 
            reads->getFlags(ad.readIds[1]).isPalindromic) continue;
        
        // Mark this alignment as in read graph
        ad.info.isInReadGraph = 1;
        
        // Create edge: (read0-strand0) to (read1-strand based on isSameStrand)
        ReadGraphEdge edge;
        edge.alignmentId = i & 0x3fff'ffff'ffff'ffff;
        edge.crossesStrands = 0;
        edge.hasInconsistentAlignment = 0;
        edge.orientedReadIds[0] = OrientedReadId(ad.readIds[0], 0);
        edge.orientedReadIds[1] = OrientedReadId(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        
        // Ensure canonical ordering
        DINARA_ASSERT(edge.orientedReadIds[0] < edge.orientedReadIds[1]);
        readGraph.edges.push_back(edge);
        
        // Also create the reverse complemented edge
        edge.orientedReadIds[0].flipStrand();
        edge.orientedReadIds[1].flipStrand();
        DINARA_ASSERT(edge.orientedReadIds[0] < edge.orientedReadIds[1]);
        readGraph.edges.push_back(edge);
    }
    
    readGraph.edges.unreserve();
    
    // Create read graph connectivity
    readGraph.connectivity.createNew(largeDataName("ReadGraphConnectivity"), largeDataPageSize);
    readGraph.connectivity.beginPass1(2 * readCount);
    for (const ReadGraphEdge& edge : readGraph.edges) {
        readGraph.connectivity.incrementCount(edge.orientedReadIds[0].getValue());
        readGraph.connectivity.incrementCount(edge.orientedReadIds[1].getValue());
    }
    readGraph.connectivity.beginPass2();
    for (size_t i = 0; i < readGraph.edges.size(); i++) {
        const ReadGraphEdge& edge = readGraph.edges[i];
        readGraph.connectivity.store(edge.orientedReadIds[0].getValue(), uint32_t(i));
        readGraph.connectivity.store(edge.orientedReadIds[1].getValue(), uint32_t(i));
    }
    readGraph.connectivity.endPass2();
    
    // Count isolated reads
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;
    for (ReadId readId = 0; readId < readCount; readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        if (readGraph.connectivity.size(orientedReadId.getValue()) == 0) {
            isolatedReadCount++;
            isolatedReadBaseCount += reads->getReadRawSequenceLength(readId);
        }
    }
    assemblerInfo->isolatedReadCount = isolatedReadCount;
    assemblerInfo->isolatedReadBaseCount = isolatedReadBaseCount;
    
    cout << timestamp << "Created read graph with " << readGraph.edges.size() << " edges." << endl;
    cout << timestamp << "Isolated reads: " << isolatedReadCount << endl;
    cout << timestamp << "createReadGraphFromFilteredAlignments completed." << endl;
}
