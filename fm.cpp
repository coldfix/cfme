// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#include <iomanip>      // setw
#include <utility>      // move

#include "number.h"
#include "fm.h"
#include "error.h"
#include "util.h"


using std::move;
using std::string;
using std::vector;
using std::endl;
using std::setw;
using std::ostream;


namespace fm
{

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
        assert_eq_size(vec.size(), num_cols);
        if (vec.empty())
            return;
        ineqs.push_back(vec.copy());
        vec.values *= -1;
        ineqs.push_back(move(vec));
    }

    void System::add_inequality(Vector&& vec)
    {
        assert_eq_size(vec.size(), num_cols);
        if (vec.empty())
            return;
        ineqs.push_back(move(vec));
    }

    Problem System::problem() const
    {
        Problem lp(num_cols);
        for (auto&& vec : ineqs) {
            lp.add_inequality(vec.values);
        }
        return lp;
    }

    // class Vector

    Vector::Vector(size_t size)
        : values(size)
    {
    }

    Vector::Vector(ValArray v)
        : values(move(v))
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
        _assert(dim >= size()<<shift, size_error);
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
        r.values = la::scaled_addition(v0.values, s0, v1.values, s1);
        return r;
    }

    ostream& operator << (ostream& o, const System& s)
    {
        for (auto&& v : s.ineqs) {
            o << v << '\n';
        }
        return o;
    }

    ostream& operator << (ostream& o, const Vector& v)
    {
        return la::print_vector(o, v.values);
    }

    bool operator == (const Vector& a, const Vector& b)
    {
        return la::equal(a.values, b.values);
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
// signigicant bits of the entropy space index. The system must have been
// created with `num_vars=nf+ni` variables.
void set_initial_state_iid(System& s, size_t nf, size_t ni)
{
    if (ni <= 1)
        return;
    size_t dim = 1<<(nf+ni);
    size_t layer1 = ((1<<ni) - 1) << nf;
    Vector v(dim);
    v.set(layer1, -1);
    for (size_t cell = 0; cell < ni; ++cell) {
        size_t var = 1 << (nf + cell);
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
// width=4, cyclic=true:
//
//     A0  A1  A2  A3
//       B0  B1  B2  B3
//
// width=4, cyclic=false:
//
//     A0  A1  A2  A3
//       B0  B1  B2
void add_causal_constraints(System& s, size_t nf, size_t ni, size_t links)
{
    size_t dim = 1<<(nf+ni);
    size_t all = dim-1;
    // for each dependent variable i, add the conditional mutual
    // independence 0 = I(i:Nd(i)|Pa(i)):
    for (size_t i = 0; i < nf; ++i) {
        size_t Var = 1<<i;
        size_t Pa = 0;
        for (size_t j = 0; j < links; ++j) {
            size_t k = (i+j) % ni;
            Pa |= 1<<(nf+k);
        }
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
        assert_eq_size(v.size(), size);
    }
    return size;
}

int get_num_vars(const Matrix& matrix)
{
    int size = get_num_cols(matrix);
    if (size == -1)
        return -1;
    _assert(is_power_of_2(size), size_error, size);
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
        lp.add_inequality(v.values);
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
        if (!lp.is_redundant(v.values))
            r.insert(r.begin() + i, move(v));
    }
    return r;
}

Vector parse_vector(string line)
{
    return la::parse_vector<Value>(line);
}

Matrix parse_matrix(const vector<string>& lines)
{
    Matrix m;
    for (auto&& v : la::parse_matrix<Value>(lines)) {
        m.push_back(move(v));
    }
    return m;
}


//----------------------------------------
// Operations
//----------------------------------------


void solve_to::run(const solve_to::Callback& cb)
{
    auto sg = cb.enter(this);
    for (int step = 0; sys.num_cols > to; ++step) {
        auto sg = cb.start_step(step);
        int best_index = to;
        int best_rank = get_rank(to);
        for (int i = to+1; i < sys.num_cols; ++i) {
            int rank = get_rank(i);
            if (rank < best_rank) {
                best_index = i;
                best_rank = rank;
            }
        }
        eliminate{sys, best_index}.run(*cb.start_eliminate(best_index));
    }
}

int solve_to::get_rank(int index) const
{
    int pos = 0;
    int neg = 0;
    for (auto&& vec : sys.ineqs) {
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

EliminatePtr solve_to::Callback::start_eliminate(int index) const
{
    return P<eliminate::Callback>(new eliminate::Callback());
}

void eliminate::run(const eliminate::Callback& cb)
{
    auto _enter = cb.enter(this);

    System s = sys.copy();

    // Partition inequality constraints into (zero, positive, negative)
    // coefficient for the given index.
    Matrix zero, pos, neg;
    for (auto&& vec : s.ineqs) {
        Value val = vec.get(index);
        if (val == 0) {
            vec.remove(index);
            zero.push_back(move(vec));
        }
        if (val > 0) {
            pos.push_back(move(vec));
        }
        if (val < 0) {
            neg.push_back(move(vec));
        }
    }

    s.ineqs = move(zero);
    --s.num_cols;

    auto _append = cb.start_append(sys.ineqs.size(), pos.size(), neg.size());

    Problem lp = s.problem();
    int i = 0;
    for (auto&& p : pos) {
        for (auto&& n : neg) {
            auto _check = cb.start_check(i++);
            Vector v = p.eliminate(n, index);
            if (!lp.is_redundant(v.values)) {
                lp.add_inequality(v.values);
                s.add_inequality(move(v));
            }
        }
    }

    sys = move(s);
}

void minimize::run(const minimize::Callback& cb)
{
    auto sg = cb.enter(this);
    fm::Problem lp = sys.problem();
    for (int i = sys.ineqs.size()-1; i >= 0; --i) {
        auto sg = cb.start_round(i);
        lp.del_row(i);
        if (lp.is_redundant(sys.ineqs[i].values)) {
            sys.ineqs.erase(sys.ineqs.begin() + i);
        }
        else {
            lp.add_inequality(sys.ineqs[i].values);
        }
    }
}


//----------------------------------------
// Status callbacks
//----------------------------------------

IO::IO(std::ostream* o, InputPtr i)
    : out(o)
    , inp(i)
{
    if (!inp) {
        inp.reset(new terminal::Input());
    }
}

// SolveTo

SG SolveToStatusOutput::enter(solve_to* ctx) const
{
    sys = &ctx->sys;
    *out << "Eliminate: " << sys->num_cols << " -> " << ctx->to << endl;
    return SG();
}

SG SolveToStatusOutput::start_step(int step) const
{
    if (inp->avail()) {
        int c = inp->get();
        if (c == 'm') {
            minimize{*sys}.run(MinimizeStatusOutput(*this));
        }
        if (c == 'c') {
            // TODO: cancel + output current progress
        }
    }
    return SG();
}

EliminatePtr SolveToStatusOutput::start_eliminate(int index) const
{
    return P<eliminate::Callback>(new EliminateStatusOutput(*this));
}

SolveToStatusOutput::~SolveToStatusOutput()
{
    *out << endl;
}

// Eliminate

SG EliminateStatusOutput::enter(eliminate* ctx) const
{
    sys = &ctx->sys;
    return SG();
}

SG EliminateStatusOutput::start_append(int z, int p, int n) const
{
    terminal::clear_current_line(*out);
    *out << "   i = " << setw(3) << sys->num_cols
        << ",  z = " << setw(4) << z
        << ",  p+n = " << setw(3) << p+n
        << "   p*n = " << setw(4) << p*n
        << std::flush;
    return SG();
}

EliminateStatusOutput::~EliminateStatusOutput()
{
    *out << endl;
}

// Minimize

SG MinimizeStatusOutput::enter(minimize* ctx) const
{
    sys = &ctx->sys;
    num_orig = ctx->sys.ineqs.size();
    return SG();
}

SG MinimizeStatusOutput::start_round(int index) const
{
    *out << "Minimizing: " << num_orig << " -> " << sys->ineqs.size()
        << "  (i=" << index << ")"
        << std::flush;
    return SG(nullptr, [this] (void*) {
        terminal::clear_current_line(*out);
    });
}

MinimizeStatusOutput::~MinimizeStatusOutput()
{
    *out << "Minimizing: " << num_orig << " -> " << sys->ineqs.size()
        << " (DONE)"
        << endl;
}


}
