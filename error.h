#ifndef __ERROR_H__INCLUDED__
#define __ERROR_H__INCLUDED__

# include "util.h"
# include <cstdlib>     // abort
# include <iostream>
# include <stdexcept>   // runtime_error


class parse_error : public std::runtime_error
{
public:
    using runtime_error::runtime_error;
};

class size_error : public std::runtime_error
{
public:
    using runtime_error::runtime_error;
};


#  define _assert(expr, exc, ...) \
    if (!(expr)) { \
        std::cerr << util::join("", \
                    __FILE__, ":", __LINE__, " in ", #expr, "\n", \
                    util::join(" ", ##__VA_ARGS__)) << std::endl; \
        abort(); \
    }

# define assert_eq(a, b, exc, ...) \
    _assert(a == b, exc, a, "!=", b, "\n", ##__VA_ARGS__);

# define assert_eq_size(a, b, ...) \
    assert_eq(a, b, size_error, ##__VA_ARGS__)

#endif // include guard
