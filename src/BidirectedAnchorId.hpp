#pragma once

// BidirectedAnchorId.hpp
//
// A bidirected anchor represents both strands of a locus as a single
// entity, matching MBG/Verkko's NodeType convention. Each bidirected
// anchor corresponds to a pair of directed anchors (anchorId, anchorId^1).
//
// BidirectedAnchorId = directed anchorId / 2.
// The strand is recovered as anchorId & 1 (0 = forward, 1 = RC).
//
// OrientedAnchor = pair<BidirectedAnchorId, bool> where bool = forward.
// This matches MBG's pair<NodeType, bool>.

#include "Shasta2Anchors.hpp"  // For Shasta2AnchorId

#include <cstdint>
#include <utility>

namespace dinara {

// Strong type to prevent mixing with directed Shasta2AnchorId.
class BidirectedAnchorId {
public:
    BidirectedAnchorId() : value_(0) {}
    explicit BidirectedAnchorId(uint64_t value) : value_(value) {}

    // Convert from a directed anchor ID.
    static BidirectedAnchorId fromDirected(Shasta2AnchorId directedId) {
        return BidirectedAnchorId(uint64_t(directedId) / 2);
    }

    // Convert back to a directed anchor ID on a given strand.
    Shasta2AnchorId toDirected(bool forward) const {
        return Shasta2AnchorId(value_ * 2 + (forward ? 0 : 1));
    }

    // The forward directed anchor (even).
    Shasta2AnchorId forwardDirected() const {
        return Shasta2AnchorId(value_ * 2);
    }

    // The RC directed anchor (odd).
    Shasta2AnchorId rcDirected() const {
        return Shasta2AnchorId(value_ * 2 + 1);
    }

    uint64_t value() const { return value_; }

    bool operator==(const BidirectedAnchorId& other) const { return value_ == other.value_; }
    bool operator!=(const BidirectedAnchorId& other) const { return value_ != other.value_; }
    bool operator<(const BidirectedAnchorId& other) const { return value_ < other.value_; }
    bool operator<=(const BidirectedAnchorId& other) const { return value_ <= other.value_; }
    bool operator>(const BidirectedAnchorId& other) const { return value_ > other.value_; }
    bool operator>=(const BidirectedAnchorId& other) const { return value_ >= other.value_; }

private:
    uint64_t value_;
};

// An oriented anchor: (BidirectedAnchorId, forward).
// Matches MBG's pair<NodeType, bool>.
using OrientedAnchor = std::pair<BidirectedAnchorId, bool>;

// Reverse an oriented anchor (flip strand).
inline OrientedAnchor reverseAnchor(OrientedAnchor a) {
    return {a.first, !a.second};
}

// Canonical form of a link between two oriented anchors.
// Matches MBG/Verkko's canon() implementation.
inline std::pair<OrientedAnchor, OrientedAnchor> canonAnchor(
    OrientedAnchor from, OrientedAnchor to)
{
    if(to.first < from.first) {
        return {reverseAnchor(to), reverseAnchor(from)};
    }
    if(to.first == from.first && !to.second && !from.second) {
        return {reverseAnchor(to), reverseAnchor(from)};
    }
    return {from, to};
}

// Check if (from, to) is the canonical direction.
inline bool isCanonicalAnchorDirection(OrientedAnchor from, OrientedAnchor to) {
    auto c = canonAnchor(from, to);
    return c.first == from && c.second == to;
}

// Convert a directed anchor ID to an oriented anchor.
inline OrientedAnchor toOrientedAnchor(Shasta2AnchorId directedId) {
    return {BidirectedAnchorId::fromDirected(directedId),
            (uint64_t(directedId) & 1) == 0};
}

} // namespace dinara

// Hash support for use in unordered containers.
namespace std {
template<> struct hash<dinara::BidirectedAnchorId> {
    size_t operator()(const dinara::BidirectedAnchorId& id) const {
        return hash<uint64_t>()(id.value());
    }
};
template<> struct hash<dinara::OrientedAnchor> {
    size_t operator()(const dinara::OrientedAnchor& oa) const {
        size_t h = hash<uint64_t>()(oa.first.value());
        h ^= hash<bool>()(oa.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
} // namespace std
