#include "HifiasmBoundaryVerify.hpp"

#include "Reads.hpp"

#include "algorithm.hpp"
#include "array.hpp"
#include "cstdint.hpp"

#include <cstring>

using namespace dinara;
using namespace std;

namespace {
    // hifiasm constants used by boundary_verify / verify_single_window.
    constexpr uint32_t HifiasmWindow = 375;
    constexpr uint32_t HifiasmThresholdMaxSize = 31;
    constexpr double HifiasmMaxOvDiffEc = 0.07;

    using Word = uint64_t;

    inline char complementBase(char c)
    {
        switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return 'N';
        }
    }

    // Recover a contiguous read subregion in either forward or reverse-complement orientation.
    inline void recoverSubregion(
        const Reads& reads,
        ReadId readId,
        uint32_t start,
        uint32_t length,
        bool reverseComplement,
        char* out)
    {
        const auto& read = reads.getRead(readId);
        const uint32_t readLen = uint32_t(read.baseCount);
        if (!reverseComplement) {
            for (uint32_t i = 0; i < length; ++i) {
                const uint32_t pos = start + i;
                out[i] = (pos < readLen) ? read[pos].character() : 'N';
            }
        } else {
            for (uint32_t i = 0; i < length; ++i) {
                const uint32_t pos = start + (length - 1U - i);
                const char b = (pos < readLen) ? read[pos].character() : 'N';
                out[i] = complementBase(b);
            }
        }
    }

    // hifiasm Adjust_Threshold macro.
    inline int adjustThreshold(int threshold, int xLen)
    {
        return (threshold == 0 && xLen >= 4) ? 1 : threshold;
    }

    // hifiasm determine_overlap_region.
    inline bool determineOverlapRegion(
        int threshold,
        int64_t yStart,
        int64_t yLen,
        int64_t windowLen,
        int& extraBegin,
        int& extraEnd,
        int64_t& yStartOut,
        int64_t& yLengthOut)
    {
        if (yStart < 0 || yLen <= yStart ||
            (yLen - yStart + 2 * int64_t(threshold) + int64_t(HifiasmThresholdMaxSize) < windowLen)) {
            return false;
        }

        extraBegin = 0;
        extraEnd = 0;
        yStartOut = yStart - threshold;
        yLengthOut = min<int64_t>(windowLen, yLen - yStartOut);
        extraEnd = int(windowLen - yLengthOut);

        if (yStartOut < 0) {
            extraBegin = int(-yStartOut);
            yStartOut = 0;
            yLengthOut -= extraBegin;
        }
        return yLengthOut > 0;
    }

    // Equivalent to fill_subregion in hifiasm.
    inline void fillSubregion(
        vector<char>& out,
        const Reads& reads,
        int64_t startPos,
        int64_t length,
        bool reverseComplement,
        ReadId readId,
        int extraBegin,
        int extraEnd)
    {
        const size_t total = size_t(extraBegin + length + extraEnd);
        out.resize(total);
        std::memset(out.data(), 'N', total);
        if (length > 0) {
            recoverSubregion(
                reads, readId, uint32_t(startPos), uint32_t(length), reverseComplement, out.data() + extraBegin);
        }
    }

    // Copy of hifiasm Reserve_Banded_BPM (Levenshtein_distance.h).
    // Returns best end-site in pattern, and sets returnErr to best error (or UINT_MAX when failing).
    inline int reserveBandedBPM(
        char* pattern,
        int pLength,
        char* text,
        int tLength,
        uint16_t errThreshold,
        unsigned int* returnErr)
    {
        (*returnErr) = (unsigned int)-1;

        Word Peq[256];

        int bandLength = (errThreshold << 1) + 1;
        int i = 0;
        Word tmpPeq = (Word)1;

        Peq[(uint8_t)'A'] = (Word)0;
        Peq[(uint8_t)'T'] = (Word)0;
        Peq[(uint8_t)'G'] = (Word)0;
        Peq[(uint8_t)'C'] = (Word)0;

        Word PeqA;
        Word PeqT;
        Word PeqC;
        Word PeqG;

        for (i = 0; i < bandLength; i++) {
            Peq[(uint8_t)pattern[i]] = Peq[(uint8_t)pattern[i]] | tmpPeq;
            tmpPeq = tmpPeq << 1;
        }

        PeqA = Peq[(uint8_t)'A'];
        PeqC = Peq[(uint8_t)'C'];
        PeqT = Peq[(uint8_t)'T'];
        PeqG = Peq[(uint8_t)'G'];

        std::memset(Peq, 0, sizeof(Word) * 256);

        Peq[(uint8_t)'A'] = PeqA;
        Peq[(uint8_t)'C'] = PeqC;
        Peq[(uint8_t)'T'] = PeqT;
        Peq[(uint8_t)'G'] = PeqG;

        Word Mask = ((Word)1 << (errThreshold << 1));

        Word VP = 0;
        Word VN = 0;
        Word X = 0;
        Word D0 = 0;
        Word HN = 0;
        Word HP = 0;

        i = 0;
        int err = 0;
        Word errMask = (Word)1;
        int iBd = (errThreshold << 1);
        int lastHigh = (errThreshold << 1);
        const int tLength1 = tLength - 1;

        while (i < tLength1) {
            X = Peq[(uint8_t)text[i]] | VN;
            D0 = ((VP + (X & VP)) ^ VP) | X;
            HN = VP & D0;
            HP = VN | ~(VP | D0);
            X = D0 >> 1;
            VN = X & HP;
            VP = HN | ~(X | HP);

            if (!(D0 & errMask)) {
                ++err;
                if ((err - lastHigh) > int(errThreshold)) {
                    return -1;
                }
            }

            Peq[(uint8_t)'A'] = Peq[(uint8_t)'A'] >> 1;
            Peq[(uint8_t)'C'] = Peq[(uint8_t)'C'] >> 1;
            Peq[(uint8_t)'G'] = Peq[(uint8_t)'G'] >> 1;
            Peq[(uint8_t)'T'] = Peq[(uint8_t)'T'] >> 1;

            ++i;
            ++iBd;
            Peq[(uint8_t)pattern[iBd]] = Peq[(uint8_t)pattern[iBd]] | Mask;
        }

        X = Peq[(uint8_t)text[i]] | VN;
        D0 = ((VP + (X & VP)) ^ VP) | X;
        HN = VP & D0;
        HP = VN | ~(VP | D0);
        X = D0 >> 1;
        VN = X & HP;
        VP = HN | ~(X | HP);
        if (!(D0 & errMask)) {
            ++err;
            if ((err - lastHigh) > int(errThreshold)) {
                return -1;
            }
        }

        int site = tLength - 1;
        int returnSite = -1;
        int availableI = pLength - tLength;
        if ((err <= int(errThreshold)) && ((unsigned int)err <= *returnErr)) {
            *returnErr = (unsigned int)err;
            returnSite = site;
        }
        i = 0;

        unsigned int ungapError = (unsigned int)-1;
        while (i < availableI) {
            err = err + int((VP >> i) & (Word)1);
            err = err - int((VN >> i) & (Word)1);
            ++i;

            if ((err <= int(errThreshold)) && ((unsigned int)err <= *returnErr)) {
                *returnErr = (unsigned int)err;
                returnSite = site + i;
            }
            if (i == int(errThreshold)) {
                ungapError = (unsigned int)err;
            }
        }

        if ((ungapError <= errThreshold) && (ungapError == (*returnErr))) {
            returnSite = site + errThreshold;
        }

        return returnSite;
    }

    // hifiasm verify_single_window equivalent.
    inline bool verifySingleWindow(
        const Reads& reads,
        uint32_t xStart,
        uint32_t xEnd,
        uint32_t overlapXs,
        uint32_t overlapYs,
        ReadId xId,
        ReadId yId,
        bool xStrand,
        vector<char>& xBuffer,
        vector<char>& yBuffer)
    {
        if (xEnd < xStart) {
            return false;
        }
        const int xLen = int(xEnd - xStart + 1);
        int threshold = int(double(xLen) * HifiasmMaxOvDiffEc);
        threshold = adjustThreshold(threshold, xLen);

        const int64_t yLen = int64_t(reads.getReadRawSequenceLength(yId));
        int64_t yStart = int64_t(xStart) - int64_t(overlapXs) + int64_t(overlapYs);
        const int64_t windowLen = int64_t(xLen) + (int64_t(threshold) << 1);

        int extraBegin = 0;
        int extraEnd = 0;
        int64_t clippedYStart = 0;
        int64_t clippedYLen = 0;
        if (!determineOverlapRegion(
                threshold, yStart, yLen, windowLen, extraBegin, extraEnd, clippedYStart, clippedYLen)) {
            return false;
        }

        // hifiasm "unusual direction":
        //   pattern is the projected y-window (always forward),
        //   text is x-subregion (in xStrand orientation).
        fillSubregion(yBuffer, reads, clippedYStart, clippedYLen, false, yId, extraBegin, extraEnd);

        xBuffer.resize(size_t(xLen));
        recoverSubregion(reads, xId, xStart, uint32_t(xLen), xStrand, xBuffer.data());

        unsigned int error = (unsigned int)-1;
        reserveBandedBPM(
            yBuffer.data(),
            int(yBuffer.size()),
            xBuffer.data(),
            int(xBuffer.size()),
            uint16_t(threshold),
            &error);

        return error != (unsigned int)-1;
    }
}

bool dinara::hifiasmBoundaryVerify(
    const Reads& reads,
    uint32_t qIntervalStart,
    uint32_t qIntervalEnd,
    ReadId qId,
    ReadId tId,
    uint32_t qs,
    uint32_t ts,
    uint32_t te,
    bool rev,
    vector<char>& xBuffer,
    vector<char>& yBuffer)
{
    // hifiasm boundary_verify mirrors:
    // dir=rev, xs=Get_qs(map), ys=(rev ? tLen-te : ts), then project query interval to target.
    if (qIntervalEnd <= qIntervalStart) {
        return false;
    }
    const uint32_t intervalLen = qIntervalEnd - qIntervalStart;

    const uint32_t tLen = uint32_t(reads.getReadRawSequenceLength(tId));
    const uint32_t xs = qs;
    const uint32_t ys = rev ? (tLen - te) : ts;

    uint32_t yIntervalStart = (qIntervalStart - xs) + ys;
    if (yIntervalStart >= tLen) {
        return false;
    }
    uint32_t yIntervalEnd = yIntervalStart + intervalLen - 1;
    if (yIntervalEnd >= tLen) {
        yIntervalEnd = tLen - 1;
    }
    if (yIntervalEnd < yIntervalStart) {
        return false;
    }

    const uint32_t yIntervalLen = yIntervalEnd - yIntervalStart + 1;
    if (yIntervalLen <= HifiasmWindow) {
        return verifySingleWindow(
            reads,
            yIntervalStart,
            yIntervalEnd,
            ys,
            xs,
            tId,
            qId,
            rev,
            xBuffer,
            yBuffer);
    }

    if (!verifySingleWindow(
            reads,
            yIntervalStart,
            yIntervalStart + HifiasmWindow - 1,
            ys,
            xs,
            tId,
            qId,
            rev,
            xBuffer,
            yBuffer)) {
        return false;
    }

    if (!verifySingleWindow(
            reads,
            yIntervalEnd - HifiasmWindow + 1,
            yIntervalEnd,
            ys,
            xs,
            tId,
            qId,
            rev,
            xBuffer,
            yBuffer)) {
        return false;
    }

    return true;
}
