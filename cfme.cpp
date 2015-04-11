#include "cfme.h"
#include "fm.h"
#include "number.h"

#include <utility>


using std::move;


// The total number of inequalities is N for the conditional entropies
// plus (N choose 2) * #(subsets) for the conditional mutual information.
template <class Int>
Int num_elemental_inequalities(Int num_vars)
{
    return num_vars + nCr(num_vars, Int(2)) * Int(1)<<(num_vars-2);
}


// Return elemental inequalities for a system of num_vars random variables.
fm::System elemental_inequalities(size_t num_vars)
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

    // Number of initial inequalities
    size_t nb_lines = num_elemental_inequalities(num_vars);

    // The first column is not used to
    //  - make the bit representation of indices match the entropy set
    //  - better integrate with GLPK's 1-based indexing
    size_t nb_cols = dim + 1;

    // Create the system
    fm::System system = fm::System(nb_lines, nb_cols);

    // index of the entropy component corresponding to the joint entropy of
    // all variables. NOTE: since the left-most column is not used, the
    // variables involved in a joint entropy correspond exactly to the bit
    // representation of its index.
    size_t all = dim;

    // Add all elemental conditional entropy positivities, i.e. those of
    // the form H(X_i|X_c)>=0 where c = ~ {i}:
    for (size_t i = 0; i < num_vars; ++i) {
        size_t c = all ^ (1 << i);
        fm::Vector v(nb_cols);
        v.set(all, 1);
        v.set(c, -1);
        system.add_inequality(move(v));
    }

    // Add all elemental conditional mutual information positivities, i.e.
    // those of the form H(X_a:X_b|X_K)>=0 where a,b not in K
    for (size_t a = 0; a < num_vars-1; ++a) {
        for (size_t b = a+1; b < num_vars; ++b) {
            size_t A = 1 << a;
            size_t B = 1 << b;
            for (size_t i = 0; i < sub_dim; ++i) {
                size_t K = skip_bit(skip_bit(i, a), b);
                fm::Vector v(nb_cols);
                v.set(A|K, 1);
                v.set(B|K, 1);
                v.set(A|B|K, -1);
                if (K) {
                    v.set(K, -1);
                }
                system.add_inequality(move(v));
            }
        }
    }

    return system;
}


// Add mutual independence constraints for the initial layer of a CCA. The
// variables of the initial layer are assumed to correspond to the most
// signigicant bits of the entropy space index. The system must be created
// with `num_vars=2*width` variables.
void set_initial_state_iid(fm::System& s, size_t width)
{
    size_t dim = 1<<(2*width);
    size_t layer1 = ((1<<width) - 1) << width;
    for (size_t cell = 0; cell < width; ++cell) {
        size_t var = 1 << (width + cell);
        size_t other = layer1 ^ var;
        fm::Vector v(dim);
        v.set(var, 1);
        v.set(other, 1);
        v.set(layer1, -1);
        s.add_equality(move(v));
    }
}


// Iterate causal constraints in first layer of a CCA of the given width.
//
// Each constraint is a conditional independency which is returned as the
// vector of its coefficients for the joint entropies.
//
// The structure of the CCA is assumed to be hexagonal:
//
//     A0  A1  A2  A3
//       B0  B1  B2  B3
void add_causal_constraints(fm::System& s, size_t width)
{
    size_t dim = 1<<(2*width);
    size_t all = dim-1;
    // for each dependent variable i, add the conditional mutual
    // independence 0 = I(i:Nd(i)|Pa(i)):
    for (size_t i = 0; i < width; ++i) {
        size_t j = (i+1) % width;
        size_t Var = 1<<i;
        size_t Pa = (1<<(width+i)) | (1<<(width+j));
        size_t Nd = all ^ (Var | Pa);
        fm::Vector v(dim);
        v.set(Pa|Var, 1);
        v.set(Pa|Nd, 1);
        v.set(Pa, -1);
        v.set(all, -1);
        s.add_equality(move(v));
    }
}


// Enumerate information inequalities in second layer of a CCA of the given
// width. The initial layer is initialized to be mutually independent. The
// layout of the CCA is as described above (c.f. `add_causal_constraints`).
bool solve(size_t width)
{
    size_t num_vars = width*2;
    size_t solve_to = 1<<width;

    fm::System system = elemental_inequalities(num_vars);
    set_initial_state_iid(system, width);
    add_causal_constraints(system, width);

    system.solve_to(solve_to);

    // used to remove inequalities implied by elemental inequalities on the
    // reduced space:
    fm::System target = elemental_inequalities(width);
    for (auto&& v : system.ineqs) {
        if (target.is_redundant(v))
            continue;
        target.add_inequality(v.copy());
        std::cout << v << std::endl;
    }
    return true;
}
