
// Helper struct for 'Flat Vector' optimization.
struct InvertedIndexTempHit {
    ReadId partnerReadId;
    uint32_t posA;
    uint32_t posB;
    
    // Sort by PartnerId then PosA.
    bool operator<(const InvertedIndexTempHit& other) const {
        if (partnerReadId != other.partnerReadId) {
            return partnerReadId < other.partnerReadId;
        }
        // Secondary sort by PosA prevents needing a second sort later?
        // Yes! If we sort list by (ID, PosA), then each group is ALREADY sorted by PosA.
        // This avoids N small sorts! Massive win.
        return posA < other.posA;
    }
};
