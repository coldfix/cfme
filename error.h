#ifndef __ERROR_H__INCLUDED__
#define __ERROR_H__INCLUDED__

# include <ostream>
# include <sstream>
# include <string>
# include <utility>

//----------------------------------------
// Compose a string message from multiple values
//----------------------------------------

inline void print_all(std::ostream& out)
{
    out << std::flush;
}

template <class H, class ...T>
void print_all(std::ostream& out, H head, T... t)
{
    out << head;
    print_all(out, std::move(t)...);
}

template <class ...T>
std::string sprint_all(T... t)
{
    std::ostringstream out;
    print_all(out, std::move(t)...);
    return out.str();
}


//----------------------------------------
// "assertion"
//----------------------------------------

template <class E, class... T>
void _assert(bool expr, T... values)
{
    if (!expr) {
        throw E(sprint_all(std::move(values)...));
    }
}

#endif // include guard
