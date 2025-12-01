#include "platformDependent.hpp"
#include <stdlib.h>
#include "fstream.hpp"

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/resource.h>
#include <mach/mach.h>
#endif

// Return the path to a usable temporary directory, including the final "/".
std::string dinara::tmpDirectory()
{
#ifdef __APPLE__
    // macOS doesn't have /dev/shm by default. /tmp is the standard.
    return "/tmp/";
#else
    // Use /dev/shm, so the temporary files never go to disk.
    return "/dev/shm/";
#endif
}



// Return the name of a timeout command or equivalent.
std::string dinara::timeoutCommand()
{
#ifdef __APPLE__
    // macOS has 'gtimeout' (from coreutils) or 'timeout' (if installed).
    // Standard macOS doesn't have a direct 'timeout' command that works like GNU timeout.
    // Assuming 'gtimeout' is available via Homebrew (coreutils), or falling back to 'perl' hack if needed.
    // For now, let's assume the user has 'gtimeout' or we return "timeout" and hope it's the GNU one.
    // A safer bet for standard macOS is to just return "timeout" and let the user install coreutils.
    return "gtimeout"; 
#else
    return "timeout";
#endif
}



uint64_t dinara::getPeakMemoryUsage() {
    uint64_t peakMemoryUsage = 0ULL;

#ifdef __APPLE__
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    // ru_maxrss is in bytes on macOS.
    peakMemoryUsage = (uint64_t)rusage.ru_maxrss;
#else
    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        while (std::getline(procStats, line)) {
            if (string::npos == line.find("VmPeak")) {
                continue;
            }
            size_t pos = line.find(":");
            while (pos < line.size() && !isdigit(line[pos])) {
                pos++;
            }
            char* end;
            peakMemoryUsage = std::strtoull(line.c_str() + pos, &end, 10);
            // Convert from kB to bytes.
            peakMemoryUsage *= 1024;
            break;
        }
    }
#endif

    return peakMemoryUsage;
}



// Get total physical memory available, in bytes.
uint64_t dinara::getTotalPhysicalMemory()
{
#ifdef __APPLE__
    int mib[2];
    int64_t physical_memory;
    size_t length;

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    length = sizeof(int64_t);
    sysctl(mib, 2, &physical_memory, &length, NULL, 0);
    return (uint64_t)physical_memory;
#else
    ifstream meminfo("/proc/meminfo");
    string s;
    uint64_t memoryKb;
    meminfo >> s >> memoryKb;
    return 1024 * memoryKb;
#endif
}

