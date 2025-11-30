// Definition of macro DINARA_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#ifndef DINARA_DINARA_ASSERT_HPP
#define DINARA_DINARA_ASSERT_HPP

namespace dinara {
    void handleFailedAssertion(
        const char* expression,
        const char* function,
        const char* file,
        int line
    ) __attribute__ ((__noreturn__));
}


#define DINARA_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (dinara::handleFailedAssertion(#expression, __PRETTY_FUNCTION__,  __FILE__ , __LINE__)))


#endif

