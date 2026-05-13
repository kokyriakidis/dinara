#ifndef DINARA_THESEUS_ALIGN_MUTEX_HPP
#define DINARA_THESEUS_ALIGN_MUTEX_HPP

#include <mutex>
#include <utility>

namespace dinara {

// theseus-lib is not thread-safe and may use process-wide state. Every
// theseus::TheseusMSA construct/align/print_as_msa sequence must hold this lock.
inline std::mutex& theseusAlignMutex()
{
    static std::mutex instance;
    return instance;
}

// Shasta2 theseusWrapper passes align(seq, weight, leftFlag, rightFlag):
//   both anchors fixed: (false, false)
//   left-anchor-only row: (false, true)
//   right-anchor-only row: (true, true)
// See https://github.com/paoloshasta/shasta2/blob/main/src/theseusWrapper.cpp
inline std::pair<bool, bool> theseusAlignEndsFreeFlags(bool hasBothAnchors, char anchorSide)
{
    if(hasBothAnchors || anchorSide == 'B') {
        return {false, false};
    }
    if(anchorSide == 'L') {
        return {false, true};
    }
    if(anchorSide == 'R') {
        return {true, true};
    }
    return {false, false};
}

} // namespace dinara

#endif
