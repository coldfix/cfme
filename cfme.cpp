// vim: path=.,/usr/include,/opt/fm-lib/include

#include <fm/system.h>
#include <fm/solution.h>
#include <fm/vector.h>

#include <fm/solver.h>

#include "cfme.h"


// binomial coefficient [n choose r]
template <class Int>
Int nCr(Int n, Int r)
{
    if (r > n)
        return 0;
    if (r * 2 > n)
        r = n - r;
    if (r == 0)
        return 1;
    Int result = n;
    for (Int i = 2; i <= r; ++i) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


template <class Int>
Int skip_bit(Int pool, size_t bit_index)
{
    Int bit = Int(1) << bit_index;
    Int left = (pool & ~(bit-1)) << 1;
    Int right = pool & (bit-1);
    return left | right;
}


s_fm_vector_t* _next_vector(s_fm_system_t* s, size_t& offs)
{
    return s->lines[offs++];
}


void _set_vector_comp(s_fm_vector_t* v, size_t i, z_type_t n)
{
    v->vector[i].num = n;
}


s_fm_system_t* elemental_inequalities(size_t num_vars)
{
    // Identify each variable with its index i from I = {0, 1, ..., N-1}.
    // Then entropy is a real valued set function from the power set of
    // indices P = 2**I. The value for the empty set can be defined to be
    // zero and is irrelevant. Therefore the dimensionality of the problem
    // is 2**N-1.
    size_t dim = (1<<num_vars) - 1;

    // After choosing 2 variables there are 2**(N-2) possible subsets of
    // the remaining N-2 variables.
    size_t sub_dim = 1 << (num_vars-2);

    // The total number of inequalities is N for the conditional entropies
    // plus (N choose 2) * #subsets for the conditional mutual information.
    size_t nb_lines = num_vars + nCr<size_t>(num_vars, 2) * sub_dim;

    // The first column signals if this is an inequality or equality. The
    // last column is the right-hand-side of the inequality.
    size_t nb_cols = dim + 2;

    // Create the system
    s_fm_system_t* system = fm_system_alloc(nb_lines, nb_cols);

    // Index of current vector.
    size_t offs = 0;

    // index of the entropy component corresponding to the joint entropy of
    // all variables. NOTE: since the left-most column is not used, the
    // variables involved in a joint entropy correspond exactly to the bit
    // representation of its index.
    size_t all = dim;

    // Add all elemental conditional entropy positivities, i.e. those of
    // the form H(X_i|X_c)>=0 where c = ~ {i}:
    for (size_t i = 0; i < num_vars; ++i) {
        size_t c = all ^ (1 << i);
        s_fm_vector_t* v = _next_vector(system, offs);
        _set_vector_comp(v, 0, 1);      // this is an inequality
        _set_vector_comp(v, all, 1);
        _set_vector_comp(v, c, -1);
        fm_vector_compute_key(&v->key, v);
    }

    // Add all elemental conditional mutual information positivities, i.e.
    // those of the form H(X_a:X_b|X_K)>=0 where a,b not in K
    for (size_t a = 0; a < num_vars-1; ++a) {
        for (size_t b = a+1; b < num_vars; ++b) {
            size_t A = 1 << a;
            size_t B = 1 << b;
            for (size_t i = 0; i < sub_dim; ++i) {
                size_t K = skip_bit(skip_bit(i, a), b);
                s_fm_vector_t* v = _next_vector(system, offs);
                _set_vector_comp(v, 0, 1);      // this is an inequality
                _set_vector_comp(v, A|K, 1);
                _set_vector_comp(v, B|K, 1);
                _set_vector_comp(v, A|B|K, -1);
                if (K) {
                    _set_vector_comp(v, K, -1);
                }
                fm_vector_compute_key(&v->key, v);
            }
        }
    }

    return system;
}


bool solve(size_t num_vars, size_t solve_to)
{
    s_fm_system_t* system = elemental_inequalities(num_vars);

    // Solve the given system. Use FM_SOLVER_AUTO_SIMPLIFY to help the solver
    // to scale automatically
    s_fm_solution_t* solution = fm_solver_solution_to(
            system, FM_SOLVER_FAST, solve_to);

    s_fm_system_t* result = fm_solution_to_system(solution);

    // Print the output solution.
    fm_system_print (stdout, result);

    // Be clean.
    fm_system_free(system);
    fm_solution_free(solution);
    fm_system_free(result);

    return true;
}
