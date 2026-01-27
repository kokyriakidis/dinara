#ifndef DINARA_ORIENTED_UNITIG_ID_HPP
#define DINARA_ORIENTED_UNITIG_ID_HPP

// Dinara.
#include "dinaraTypes.hpp"

// Standard library.
#include "cstdint.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

namespace dinara {
    class OrientedUnitigId;
}


class dinara::OrientedUnitigId {
public:
    OrientedUnitigId() = default;

    OrientedUnitigId(uint32_t unitigId, Strand strand) :
        unitigId_(unitigId),
        strand_(strand)
    {}

    explicit OrientedUnitigId(const std::string& s)
    {
        const auto dash = s.find('-');
        if (dash == std::string::npos) {
            throw std::runtime_error("Invalid oriented unitig id.");
        }
        unitigId_ = uint32_t(std::stoul(s.substr(0, dash)));
        strand_ = Strand(std::stoul(s.substr(dash + 1)));
        if (strand_ != 0 && strand_ != 1) {
            throw std::runtime_error("Invalid oriented unitig strand.");
        }
    }

    uint32_t getUnitigId() const { return unitigId_; }
    Strand getStrand() const { return strand_; }
    uint32_t getValue() const { return (unitigId_ << 1) | uint32_t(strand_); }
    std::string getString() const { return std::to_string(unitigId_) + "-" + std::to_string(uint32_t(strand_)); }

    static OrientedUnitigId fromValue(uint32_t orientedValue)
    {
        return OrientedUnitigId(orientedValue >> 1, Strand(orientedValue & 1U));
    }

    friend bool operator<(const OrientedUnitigId& a, const OrientedUnitigId& b)
    {
        if (a.unitigId_ != b.unitigId_) return a.unitigId_ < b.unitigId_;
        return a.strand_ < b.strand_;
    }
    friend bool operator==(const OrientedUnitigId& a, const OrientedUnitigId& b)
    {
        return a.unitigId_ == b.unitigId_ && a.strand_ == b.strand_;
    }

private:
    uint32_t unitigId_ = 0;
    Strand strand_ = 0;
};

#endif
