// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

#include "performanceLog.hpp"

namespace dinara {
    ofstream performanceLog;
}



void dinara::openPerformanceLog(const string& fileName)
{
    performanceLog.open(fileName);
}
