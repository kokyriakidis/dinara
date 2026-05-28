#ifndef DINARA_SDUST_HPP
#define DINARA_SDUST_HPP

// Symmetric DUST low-complexity filter.
// C++ port of Heng Li's sdust.c (minimap2).
// Detects low-complexity regions in 2-bit encoded DNA
// (A=0, C=1, G=2, T=3, >=4 = N/break).
// Returns half-open intervals [start, end).

#include <cstdint>
#include <cstring>
#include <deque>
#include <utility>
#include <vector>

namespace dinara {

namespace sdust_detail {

static constexpr int WLEN = 3;
static constexpr int WTOT = (1 << (WLEN << 1));  // 64
static constexpr int WMSK = WTOT - 1;

struct PerfIntv { int start, finish, r, l; };

inline void shiftWindow(
    int t, std::deque<int>& w, int T, int W,
    int& L, int& rw, int& rv, int* cw, int* cv)
{
    int s;
    if(int(w.size()) >= W - WLEN + 1) {
        s = w.front(); w.pop_front();
        rw -= --cw[s];
        if(L > int(w.size()))
            --L, rv -= --cv[s];
    }
    w.push_back(t);
    ++L;
    rw += cw[t]++;
    rv += cv[t]++;
    if(cv[t] * 10 > T << 1) {
        do {
            s = w[w.size() - L];
            rv -= --cv[s];
            --L;
        } while(s != t);
    }
}

inline void saveMasked(
    std::vector<std::pair<uint32_t,uint32_t>>& res,
    std::vector<PerfIntv>& P, int start)
{
    if(P.empty() || P.back().start >= start) return;
    auto& p = P.back();
    int saved = 0;
    if(!res.empty()) {
        uint32_t f = res.back().second;
        if(uint32_t(p.start) <= f) {
            saved = 1;
            res.back().second =
                (f > uint32_t(p.finish)) ? f : uint32_t(p.finish);
        }
    }
    if(!saved)
        res.push_back({uint32_t(p.start), uint32_t(p.finish)});
    int i;
    for(i = int(P.size()) - 1;
        i >= 0 && P[i].start < start; --i);
    P.resize(i + 1);
}

inline void findPerfect(
    std::vector<PerfIntv>& P, const std::deque<int>& w,
    int T, int start, int L, int rv, const int* cv)
{
    int c[WTOT], r = rv;
    std::memcpy(c, cv, WTOT * sizeof(int));
    int max_r = 0, max_l = 0;
    for(int i = int(w.size()) - L - 1; i >= 0; --i) {
        int t = w[i];
        r += c[t]++;
        int new_r = r, new_l = int(w.size()) - i - 1;
        if(new_r * 10 > T * new_l) {
            int j;
            for(j = 0; j < int(P.size())
                && P[j].start >= i + start; ++j) {
                auto& p = P[j];
                if(max_r == 0 || p.r * max_l > max_r * p.l)
                    max_r = p.r, max_l = p.l;
            }
            if(max_r == 0
               || new_r * max_l >= max_r * new_l) {
                max_r = new_r; max_l = new_l;
                PerfIntv pi;
                pi.start = i + start;
                pi.finish = int(w.size()) + (WLEN - 1) + start;
                pi.r = new_r; pi.l = new_l;
                P.insert(P.begin() + j, pi);
            }
        }
    }
}

} // namespace sdust_detail

// Run SDUST on a 2-bit encoded sequence.
// seq[i] must be 0-3 for valid bases, >=4 for N/break.
// T = complexity threshold (default 20, lower = more aggressive masking).
// W = window size (default 64).
// Returns sorted, non-overlapping half-open intervals.
inline void sdust(
    const uint8_t* seq, uint32_t len,
    int T, int W,
    std::vector<std::pair<uint32_t,uint32_t>>& out)
{
    using namespace sdust_detail;
    out.clear();
    if(int(len) < WLEN) return;

    std::deque<int> w;
    std::vector<PerfIntv> P;
    int rv = 0, rw = 0, L = 0;
    int cv[WTOT], cw[WTOT];
    std::memset(cv, 0, sizeof(cv));
    std::memset(cw, 0, sizeof(cw));
    unsigned t = 0;
    int l = 0;

    for(uint32_t i = 0; i <= len; ++i) {
        int b = (i < len) ? int(seq[i]) : 4;
        if(b < 4) {
            ++l;
            t = (t << 2 | b) & WMSK;
            if(l >= WLEN) {
                int start = (l - W > 0 ? l - W : 0)
                            + int(i + 1) - l;
                saveMasked(out, P, start);
                shiftWindow(t, w, T, W, L, rw, rv, cw, cv);
                if(rw * 10 > L * T)
                    findPerfect(P, w, T, start, L, rv, cv);
            }
        } else {
            int start = (l - W + 1 > 0 ? l - W + 1 : 0)
                        + int(i + 1) - l;
            while(!P.empty()) saveMasked(out, P, start++);
            l = 0; t = 0;
        }
    }
}

// Convenience overload for Base sequences (Base.value is 0-3).
template<typename BaseType>
inline void sdust(
    const std::vector<BaseType>& bases,
    uint32_t start, uint32_t end,
    int T, int W,
    std::vector<std::pair<uint32_t,uint32_t>>& out)
{
    if(end <= start) { out.clear(); return; }
    const uint32_t len = end - start;
    std::vector<uint8_t> enc(len);
    for(uint32_t i = 0; i < len; ++i)
        enc[i] = bases[start + i].value;
    sdust(enc.data(), len, T, W, out);
    // Shift intervals back to original coordinates.
    for(auto& [s, e] : out) {
        s += start;
        e += start;
    }
}

} // namespace dinara

#endif // DINARA_SDUST_HPP
