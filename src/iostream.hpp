#ifndef DINARA_IOSTREAM_HPP
#define DINARA_IOSTREAM_HPP

#include <iostream>

namespace dinara {
    using std::cin;
    using std::cout;
    using std::dec;
    using std::endl;
    using std::flush;
    using std::hex;
    using std::istream;
    using std::ostream;

    // In Dinara we don't use cerr. All log output is to cout.
    // using std::cerr;

    // Output 128-bit integer.
    inline ostream& operator<<(ostream& s, __uint128_t)
    {
        s << "(__uint128_t)";
        return s;
    }
}

#endif
