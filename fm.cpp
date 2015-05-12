// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).


#include <algorithm>    // copy
#include <cassert>
#include <cmath>
#include <iomanip>      // setw
#include <iterator>     // istream_iterator, back_inserter
#include <utility>      // move

#include <glpk.h>

#include "number.h"
#include "fm.h"
#include "error.h"
#include "util.h"


using std::move;
using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::setw;


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

    void Problem::del_row(int i)
    {
        glp_del_rows(prob.get(), 1, (&i)-1);
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
            throw std::runtime_error("Error in glp_simplex.");
        }
        return glp_get_status(prob.get()) == GLP_OPT;
    }

    bool Problem::dual(const Vector& v, std::vector<double>& r) const
    {
        assert(v.size() == num_cols);
        glp_prob* lp = prob.get();
        for (int i = 1; i < num_cols; ++i) {
            glp_set_obj_coef(lp, i, v.get(i));
        }
        glp_std_basis(lp);
        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_ERR;
        parm.meth = GLP_DUAL;
        int result = glp_simplex(lp, &parm);
        if (result != 0) {
            throw std::runtime_error("Error in glp_simplex.");
        }
        int status = glp_get_dual_stat(lp);
        if (status != GLP_FEAS) {
            return false;
        }
        int rows = glp_get_num_cols(lp);
        r.clear();
        r.resize(rows);
        for (int i = 0; i < rows; ++i) {
            r[i] = glp_get_row_dual(lp, i+1);
        }
        return true;
    }

    // class System

    System::System(size_t nb_lines, size_t nb_cols)
        : num_cols(nb_cols)
    {
        clear(nb_lines);
    }

    System::System(Matrix matrix)
    {
        num_cols = get_num_cols(matrix);
        ineqs = move(matrix);
    }

    System System::copy() const
    {
        System s(ineqs.size(), num_cols);
        for (auto&& v : ineqs) {
            s.add_inequality(v.copy());
        }
        return s;
    }

    void System::clear(size_t new_expected)
    {
        ineqs.clear();
        ineqs.reserve(new_expected);
    }

    void System::add_equality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (vec.empty())
            return;
        ineqs.push_back(vec.copy());
        vec.values *= -1;
        ineqs.push_back(move(vec));
    }

    void System::add_inequality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (vec.empty())
            return;
        ineqs.push_back(move(vec));
    }

    void System::solve_to(int to, int* recorded_order)
    {
        int num_orig = num_cols;
        cerr << "Eliminate: " << num_orig << " -> " << to << endl;
        for (int step = 0; num_cols > to; ++step) {
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
            if (recorded_order) {
                recorded_order[step] = best_index;
            }
        }
        cerr << endl;
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
        Problem lp(num_cols);
        for (auto&& vec : ineqs) {
            lp.add_inequality(vec);
        }
        return lp;
    }

    void System::eliminate(int index)
    {
        int num_orig_ineqs = ineqs.size();

        --num_cols;
        Matrix _ineqs = move(ineqs);
        clear(_ineqs.size());

        // Partition inequality constraints into (positive, negative, zero)
        // coefficient for the given index.
        Matrix pos, neg;
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

        Matrix cand;
        cand.reserve(pos.size()*neg.size());

        for (auto&& p : pos) {
            for (auto&& n : neg) {
                cand.push_back(p.eliminate(n, index));
            }
        }

        Problem lp = problem();

        terminal::clear_current_line(cerr);
        cerr
            << "   i = " << setw(3) << num_cols
            << ",  num_ineqs = " << setw(4) << ineqs.size()
            << ",  p+n = " << setw(3) << pos.size()+neg.size()
            << "   p*n = " << setw(4) << pos.size()*neg.size()
            << std::flush;
        for (int i = 0; i < cand.size(); ++i) {
            Vector& vec = cand[i];
            if (!lp.is_redundant(vec)) {
                lp.add_inequality(vec);
                add_inequality(move(vec));
            }
        }
        cerr << endl;
        if (ineqs.size() > num_orig_ineqs + 10) {
            minimize();
            terminal::cursor_up(cerr);
            terminal::clear_current_line(cerr);
        }
    }

    void System::minimize()
    {
        size_t num_orig = ineqs.size();
        fm::Problem lp = problem();
        for (int i = ineqs.size()-1; i >= 0; --i) {
            cerr << "Minimizing: " << num_orig << " -> " << ineqs.size()
                << "  (i=" << i << ")"
                << std::flush;
            lp.del_row(i+1);
            if (lp.is_redundant(ineqs[i])) {
                ineqs.erase(ineqs.begin() + i);
            }
            else {
                lp.add_inequality(ineqs[i]);
            }
            terminal::clear_current_line(cerr);
        }
        cerr << "Minimizing: " << num_orig << " -> " << ineqs.size()
            << " (DONE)"
            << endl;
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

    // friends & co

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
            o << setw(3) << val << ' ';
        }
        o << "]";
        return o;
    }

    bool operator == (const Vector& a, const Vector& b)
    {
        assert(a.size() == b.size());
        for (int i = 0; i < a.size(); ++i) {
            if (a.get(i) != b.get(i)) {
                return false;
            }
        }
        return true;
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
    if (width <= 1)
        return;
    size_t dim = 1<<(2*width);
    size_t layer1 = ((1<<width) - 1) << width;
    Vector v(dim);
    v.set(layer1, -1);
    for (size_t cell = 0; cell < width; ++cell) {
        size_t var = 1 << (width + cell);
        v.set(var, 1);
    }
    s.add_equality(move(v));
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


//----------------------------------------
// matrix "methods"
//----------------------------------------

int get_num_cols(const Matrix& matrix)
{
    if (matrix.empty())
        return -1;
    int size = matrix[0].size();
    for (auto&& v : matrix) {
        _assert<matrix_size_error>(v.size() == size,
                "size does not match", v.copy());
    }
    return size;
}

int get_num_vars(const Matrix& matrix)
{
    int size = get_num_cols(matrix);
    if (size == -1)
        return -1;
    _assert<matrix_size_error>(is_power_of_2(size),
            "size must be power of 2", size);
    return intlog2(size);
}

Matrix copy_matrix(const Matrix& m)
{
    Matrix r;
    for (auto&& v : m)
        r.push_back(v.copy());
    return r;
}

Problem problem(const Matrix& m, int num_vars)
{
    fm::System sys(m.size(), 1<<num_vars);
    fm::Problem lp = sys.problem();
    for (auto&& v : m)
        lp.add_inequality(v.copy());
    return lp;
}

// greedy minimization
Matrix minimize_system(const Matrix& sys)
{
    int num_vars = get_num_vars(sys);
    Matrix r = copy_matrix(sys);
    for (int i = r.size()-1; i > 0; --i) {
        Vector v = move(r[i]);
        r.erase(r.begin() + i);
        Problem lp = problem(r, num_vars);
        if (!lp.is_redundant(v))
            r.insert(r.begin() + i, move(v));
    }
    return r;
}

string trim(string s)
{
    int beg = s.find_first_not_of(" \t");
    int end = s.find_last_not_of(" \t");
    if (beg == -1)
        return string();
    return s.substr(beg, end-beg+1);
}

string remove_comment(string s)
{
    int beg = s.find('#');
    if (beg == -1)
        return s;
    return s.substr(0, beg);
}

Vector parse_vector(string line)
{
    typedef std::istream_iterator<int> iit;
    _assert<matrix_parse_error>(line.front() == '[', "expecting '['", line);
    _assert<matrix_parse_error>(line.back() == ']', "expecting ']'", line);
    line = trim(line.substr(1, line.size()-2));
    std::istringstream in(line);
    vector<int> vals;
    copy(iit(in), iit(), std::back_inserter(vals));
    Vector r(vals.size());
    copy(vals.begin(), vals.end(), begin(r.values));
    return r;
}

Matrix parse_matrix(const vector<string>& lines)
{
    Matrix r;
    for (string line : lines) {
        line = remove_comment(line);
        line = trim(line);
        if (line.empty())
            continue;
        r.push_back(parse_vector(line));
    }
    get_num_cols(r);
    return r;
}


}
