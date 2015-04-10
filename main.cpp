#include <cstdlib>          // atol
#include "cfme.h"

int main(int argc, char** argv, char** env)
{
    size_t num_vars = 2;
    size_t solve_to = 2;
    if (argc >= 2)
        num_vars = std::atol(argv[1]);
    if (argc >= 3)
        solve_to = std::atol(argv[2]);
    solve(num_vars, solve_to);
    return 0;
}
