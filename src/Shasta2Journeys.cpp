#include "Shasta2Journeys.hpp"
#include "Shasta2Anchors.hpp"
#include "DINARA_ASSERT.hpp"
#include "orderPairs.hpp"
#include "MultithreadedObject.tpp"
#include "Reads.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace dinara;
using namespace std;

Shasta2Journeys::Shasta2Journeys(
    uint64_t orientedReadCount,
    shared_ptr<Shasta2Anchors> anchorsPointer,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    anchorsPointer(anchorsPointer)
{
    // Pass 1: count occurrences to determine sizes for temporary storage.
    // Use journeysWithOrdinals to store (Shasta2AnchorId, ordinal) for each oriented read.
    // We need to size it properly.
    // VectorOfVectors needs 2 passes: 1 to count, 2 to fill.
    // journeysWithOrdinals is VectorOfVectors<pair<uint64_t, uint32_t>, uint64_t>
    
    journeysWithOrdinals.createNew(
        largeDataName("Shasta2Journeys-WithOrdinals"), 
        largeDataPageSize);
        
    // Pass 1: count how many anchors each oriented read visits.
    journeysWithOrdinals.beginPass1(orientedReadCount);
    
    const uint64_t anchorCount = anchorsPointer->size();
    const uint64_t batchSize = 1000;
    
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&Shasta2Journeys::threadFunction1, threadCount);
    
    // Pass 2: store (Shasta2AnchorId, ordinal)
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&Shasta2Journeys::threadFunction2, threadCount);
    journeysWithOrdinals.endPass2();
    
    // Pass 3: Sort each vector in journeysWithOrdinals and size 'journeys'.
    journeys.createNew(
        largeDataName("Shasta2Journeys-Journeys"),
        largeDataPageSize);
        
    journeys.beginPass1(orientedReadCount);
    setupLoadBalancing(orientedReadCount, batchSize);
    runThreads(&Shasta2Journeys::threadFunction3, threadCount); // This sorts and counts for next stage.
    
    // Pass 4: Copy sorted AnchorIds to 'journeys' and update Anchors.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, batchSize); 
    runThreads(&Shasta2Journeys::threadFunction4, threadCount);
    journeys.endPass2(false, true); // (optimize, hugepages)
    
    journeysWithOrdinals.remove();
}

Shasta2Journeys::Shasta2Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Shasta2Journeys-Journeys"));
}


void Shasta2Journeys::threadFunction1(uint64_t threadId) {
    threadFunction12(1);
}

void Shasta2Journeys::threadFunction2(uint64_t threadId) {
    threadFunction12(2);
}

void Shasta2Journeys::threadFunction12(uint64_t pass) {
    Shasta2Anchors& anchors = *anchorsPointer; // non-const access? 
    // Actually we only read from anchors here.
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(Shasta2AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            // Iterate over markers in this anchor.
            // Accessing anchors[anchorId] gives Shasta2Anchor which is span<const Shasta2AnchorMarkerInfo>.
            // We need to know orientedReadId and ordinal.
            
            const auto anchor = anchors[anchorId];
            for(const auto& info : anchor) {
                uint64_t orientedReadIdValue = info.orientedReadId.getValue();
                
                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, 
                        make_pair(uint64_t(anchorId), info.ordinal));
                }
            }
        }
    }
}

void Shasta2Journeys::threadFunction3(uint64_t threadId) {
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t orVal=begin; orVal!=end; orVal++) {
            // Sort the vector for this oriented read.
            // journeysWithOrdinals gives access to a mutable span?
            // VectorOfVectors::operator[] returns span.
            // But we need to toggle write access? it is open for write/modify?
            // "journeysWithOrdinals" is just created.
            // Wait, VectorOfVectors only gives const span via operator[] usually?
            // Need to check if we can modify elements in place.
            // Shasta2 uses `auto v = journeysWithOrdinals[orientedReadValue];` and calls `sort`.
            // So `operator[]` returns a span that allows modification if the underlying data is mutable.
            // Since we created it, it should be mutable.
            
            auto v = journeysWithOrdinals[orVal];
            
            // Sort by ordinal (second element of pair).
            // Using lambda if orderPairs.hpp not usable directly or different types to Shasta2.
            std::sort(v.begin(), v.end(), [](const pair<uint64_t, uint32_t>& a, const pair<uint64_t, uint32_t>& b) {
                return a.second < b.second;
            });
            
            // Count for journeys (size matches).
            journeys.incrementCountMultithreaded(orVal, v.size());
        }
    }
}

void Shasta2Journeys::threadFunction4(uint64_t threadId) {
    Shasta2Anchors& anchors = *anchorsPointer;
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t orVal=begin; orVal!=end; orVal++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orVal));
            
            const auto src = journeysWithOrdinals[orVal]; // sorted
            auto dest = journeys[orVal]; // span for writing (journeys is in pass 2)
            
            // Copy AnchorIds (first of pair) to dest.
            for(size_t i=0; i<src.size(); i++) {
                dest[i] = Shasta2AnchorId(src[i].first); // Cast to Shasta2AnchorId
            }
            
            // Update positionInJourney in Anchors.
            // We iterate through the journey we just built.
            for(size_t position=0; position<dest.size(); position++) {
                Shasta2AnchorId anchorId = dest[position];
                
                // We need to modify the Shasta2AnchorMarkerInfo in anchors.
                // anchors[anchorId] returns Shasta2Anchor (const span).
                // We need mutable access to anchors.anchorMarkerInfos[anchorId].
                // anchors.anchorMarkerInfos is accessible (public in Shasta2Anchors).
                
                auto markerInfos = anchors.anchorMarkerInfos[anchorId]; // returns mutable span? 
                // Checks MemoryMappedVectorOfVectors.hpp for non-const operator[].
                // Assuming yes.
                
                bool found = false;
                for(auto& markerInfo : markerInfos) {
                    if(markerInfo.orientedReadId == orientedReadId) {
                        markerInfo.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                // Assert found?
                DINARA_ASSERT(found);
            }
        }
    }
}
