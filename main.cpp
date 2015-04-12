#include <cstdlib>          // atol
#include "cfme.h"

int main(int argc, char** argv, char** env)
try
{
    size_t width = 2;
    if (argc >= 2)
        width = std::atol(argv[1]);
    solve(width);
    return 0;
}
catch (...)
{
    throw;
}
