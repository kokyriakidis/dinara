#include <stdexcept>
#include <string>

namespace shasta2 {

void handleFailedAssertion(
    const char* expression,
    const char* function,
    const char* file,
    int line)
{
    throw std::runtime_error(
        std::string("SHASTA2_ASSERT failed: ") + expression +
        " in " + function +
        " at " + file + ":" + std::to_string(line));
}

}
