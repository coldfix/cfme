#ifndef __ERROR_H__INCLUDED__
#define __ERROR_H__INCLUDED__

# include <utility>     // move
# include "util.h"

template <class E, class... T>
void _assert(bool expr, T... values)
{
    if (!expr) {
        throw E(util::sprint_all(std::move(values)...));
    }
}

#endif // include guard
