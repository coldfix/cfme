// Thin C++ wrapper classes around the API of the FM library.


#include <cassert>
#include <cmath>
#include <utility>      // move

#include <glpk.h>

#include "number.h"
#include "fm.h"


using std::move;


namespace fm
{

    // class Problem

    Problem::Problem()
    {
    }

    Problem::Problem(size_t nb_cols)
        : num_cols(nb_cols)
    {
        prob.reset(glp_create_prob(), glp_delete_prob);
        glp_set_obj_dir(prob.get(), GLP_MIN);
        glp_add_cols(prob.get(), num_cols-1);
        for (int j = 1; j < num_cols; ++j) {
            glp_set_col_bnds(prob.get(), j, GLP_FR, NAN, NAN);
        }
    }

    void Problem::set_mat_row(int i, const Vector& v)
    {
        std::vector<int> indices;
        std::vector<double> values;
        indices.reserve(v.size());
        values.reserve(v.size());
        indices.push_back(0);       // ind[0] is not used by GLPK
        values.push_back(NAN);      // val[0] is not used by GLPK
        for (int i = 1; i < v.size(); ++i) {
            Value val = v.get(i);
            if (val) {
                indices.push_back(i);
                values.push_back(val);
            }
        }
        glp_set_mat_row(prob.get(), i,
                indices.size()-1, indices.data(), values.data());
    }

    void Problem::add_equality(const Vector& v)
    {
        int i = glp_add_rows(prob.get(), 1);
        glp_set_row_bnds(prob.get(), i, GLP_FX, 0.0, 0.0);
        set_mat_row(i, v);
    }

    void Problem::add_inequality(const Vector& v)
    {
        int i = glp_add_rows(prob.get(), 1);
        glp_set_row_bnds(prob.get(), i, GLP_LO, 0.0, NAN);
        set_mat_row(i, v);
    }

    bool Problem::is_redundant(const Vector& v) const
    {
        assert(v.size() == num_cols);
        for (int i = 1; i < num_cols; ++i) {
            glp_set_obj_coef(prob.get(), i, v.get(i));
        }
        glp_std_basis(prob.get());
        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_ERR;
        int result = glp_simplex(prob.get(), &parm);
        if (result != 0) {
            return false;       // TODO: ERROR, raise exception?
        }
        return glp_get_status(prob.get()) == GLP_OPT;
    }

    // class System

    System::System(size_t nb_lines, size_t nb_cols)
        : num_cols(nb_cols)
    {
        clear(nb_lines);
    }

    System System::copy() const
    {
        System s(ineqs.size(), num_cols);
        for (auto&& v : ineqs) {
            s.add_inequality(v.copy());
        }
        for (auto&& v : eqns) {
            s.add_equality(v.copy());
        }
        return s;
    }

    void System::clear(size_t new_expected)
    {
        ineqs.clear();
        eqns.clear();
        ineqs.reserve(new_expected);
    }

    void System::add_equality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (vec.empty())
            return;
        eqns.push_back(move(vec));
    }

    void System::add_inequality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (vec.empty())
            return;
        ineqs.push_back(move(vec));
    }

    void System::solve_to(int to)
    {
        while (num_cols > to) {
            int best_index = to;
            int best_rank = get_rank(to);
            for (int i = to+1; i < num_cols; ++i) {
                int rank = get_rank(i);
                if (rank < best_rank) {
                    best_index = i;
                    best_rank = rank;
                }
            }
            eliminate(best_index);
        }

    }

    int System::get_rank(int index) const
    {
        int pos = 0;
        int neg = 0;
        for (auto&& vec : ineqs) {
            Value val = vec.get(index);
            if (val > 0) {
                ++pos;
            }
            else if (val < 0) {
                ++neg;
            }
        }
        return (pos*neg) - (pos+neg);
    }

    Problem System::problem() const
    {
        Problem p(num_cols);
        for (auto&& vec : eqns) {
            p.add_equality(vec);
        }
        for (auto&& vec : ineqs) {
            p.add_inequality(vec);
        }
        return p;
    }

    void System::eliminate(int index)
    {
        --num_cols;
        std::vector<Vector> _ineqs = move(ineqs);
        std::vector<Vector> _eqns = move(eqns);
        clear(_ineqs.size());

        // Partition inequality constraints into (positive, negative, zero)
        // coefficient for the given index.
        std::vector<Vector> pos, neg;
        for (auto&& vec : _ineqs) {
            Value val = vec.get(index);
            if (val > 0) {
                pos.push_back(move(vec));
            }
            else if (val < 0) {
                neg.push_back(move(vec));
            }
            else if (val == 0) {
                vec.remove(index);
                add_inequality(move(vec));
            }
        }

        std::cout
            << "cols: " << num_cols
            << ", ineqs: " << ineqs.size()
            << ", eqs: " << _eqns.size()
            << ", p*n: " << pos.size()*neg.size()
            << ", p+n: " << pos.size()+neg.size()
            << std::endl;

        std::vector<Vector> eq_with;
        eq_with.reserve(_eqns.size());
        for (auto&& vec : _eqns) {
            if (vec.get(index)) {
                eq_with.push_back(move(vec));
            }
            else {
                vec.remove(index);
                add_equality(move(vec));
            }
        }

        if (!eq_with.empty()) {
            // TODO: heuristic for choosing equation?
            Vector eq = move(eq_with.back());
            eq_with.pop_back();
            for (auto&& v : eq_with) {
                add_equality(v.eliminate(eq, index));
            }
            for (auto&& v : pos) {
                add_inequality(v.eliminate(eq, index));
            }
            for (auto&& v : neg) {
                add_inequality(v.eliminate(eq, index));
            }
        }
        else {
            Problem prob = problem();

            for (auto&& p : pos) {
                for (auto&& n : neg) {
                    Vector v = p.eliminate(n, index);
                    if (!prob.is_redundant(v)) {
                        prob.add_inequality(v);
                        add_inequality(move(v));
                    }
                }
            }
        }
    }


    // class Vector

    Vector::Vector(size_t size)
        : values(size)
    {
    }

    Vector Vector::copy() const
    {
        return Vector(*this);
    }

    bool Vector::empty() const
    {
        for (const Value& x : values) {
            if (x) {
                return false;
            }
        }
        return true;
    }

    size_t Vector::size() const
    {
        return values.size();
    }

    Value Vector::get(size_t i) const
    {
        return values[i];
    }

    void Vector::set(size_t i, Value v)
    {
        values[i] = v;
    }

    void Vector::remove(size_t i)
    {
        size_t s = size();
        std::valarray<Value> r(s - 1);
        for (int j = 0; j < i; ++j) {
            r[j] = values[j];
        }
        for (int j = i+1; j < s; ++j) {
            r[j-1] = values[j];
        }
        values.swap(r);
    }

    // Inplace elimination of coefficients.
    Vector Vector::eliminate(const Vector& v, size_t i) const
    {
        Value a = get(i);
        Value b = v.get(i);
        Value s = -sign(a*b);
        a = abs(a);
        b = abs(b);
        Value div = gcd(a, b);
        Vector r = scaled_addition(*this, (b / div),
                                   v, s * (a / div));
        r.normalize();
        r.remove(i);
        return r;
    }

    // Inplace normalization of coefficients.
    void Vector::normalize()
    {
        Value div(0);
        for (const Value& x : values) {
            div = gcd<Value>(div, abs(x));
            if (div == 1) {
                return;
            }
        }
        if (div > 1) {
            values /= div;
        }
    }

    Vector Vector::injection(size_t dim, size_t shift) const
    {
        assert(dim >= size()<<shift);
        Vector r(dim);
        for (size_t i = 0; i < size(); ++i) {
            r.set(i<<shift, get(i));
        }
        return r;
    }

    // friends

    Vector scaled_addition(const Vector& v0, Value s0,
                           const Vector& v1, Value s1)
    {
        Vector r;
        r.values = scaled_addition(v0.values, s0, v1.values, s1);
        return r;
    }

    ValArray scaled_addition(const ValArray& v0, Value s0,
                             const ValArray& v1, Value s1)
    {
        assert(v0.size() == v1.size());
        return v0 * s0 + v1 * s1;
    }

    std::ostream& operator << (std::ostream& o, const System& s)
    {
        for (auto&& v : s.ineqs) {
            o << v << '\n';
        }
        return o;
    }

    std::ostream& operator << (std::ostream& o, const Vector& v)
    {
        o << "[ ";
        for (auto val : v.values) {
            o << val << ' ';
        }
        o << "]";
        return o;
    }


//----------------------------------------
// stand-alone functions
//----------------------------------------

// The total number of inequalities is N for the conditional entropies
// plus (N choose 2) * #(subsets) for the conditional mutual information.
size_t num_elemental_inequalities(size_t num_vars)
{
    return num_vars + nCr(num_vars, 2ul) * (1ul<<(num_vars-2));
}


// Return elemental inequalities for a system of num_vars random variables.
System elemental_inequalities(size_t num_vars)
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
    System system = System(nb_lines, nb_cols);

    // index of the entropy component corresponding to the joint entropy of
    // all variables. NOTE: since the left-most column is not used, the
    // variables involved in a joint entropy correspond exactly to the bit
    // representation of its index.
    size_t all = dim;

    // Add all elemental conditional entropy positivities, i.e. those of
    // the form H(X_i|X_c)>=0 where c = ~ {i}:
    for (size_t i = 0; i < num_vars; ++i) {
        size_t c = all ^ (1 << i);
        Vector v(nb_cols);
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
                Vector v(nb_cols);
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
void set_initial_state_iid(System& s, size_t width)
{
    size_t dim = 1<<(2*width);
    size_t layer1 = ((1<<width) - 1) << width;
    for (size_t cell = 0; cell < width; ++cell) {
        size_t var = 1 << (width + cell);
        size_t other = layer1 ^ var;
        Vector v(dim);
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
void add_causal_constraints(System& s, size_t width)
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
        Vector v(dim);
        v.set(Pa|Var, 1);
        v.set(Pa|Nd, 1);
        v.set(Pa, -1);
        v.set(all, -1);
        s.add_equality(move(v));
    }
}


}
