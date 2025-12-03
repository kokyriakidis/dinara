

#include <cstdint>
#include "touchMemory.hpp"
using namespace dinara;

// Touch a range of memory in order to cause the
// supporting pages of virtual memory to be loaded in real memory.
// The return value can be ignored.
size_t dinara::touchMemory(
    const void* begin,
    const void* end,
    size_t pageSize)
{
    const char* cBegin = static_cast<const char*>(begin);
    const char* cEnd = static_cast<const char*>(end);
    size_t sum = 0ULL;
    for(const char* p=cBegin; p<cEnd; p+=pageSize) {
        sum += static_cast<uint8_t>(*p);
    }
    return sum;
}

