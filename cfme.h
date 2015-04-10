// Enumerate information inequalities in subsequent layers of a periodic CCA
// with the given size.

#ifndef __CFME_H__INCLUDED__
#define __CFME_H__INCLUDED__

# include <cstddef>

// num_vars: total number of cells (=2*width)
// solve_to: number of vars after elimination
bool solve(size_t num_vars, size_t solve_to);

#endif  // include guard
