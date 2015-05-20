// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#ifndef __FM_H__INCLUDED__
#define __FM_H__INCLUDED__

# include <iostream>
# include <memory>
# include <stdexcept>
# include <valarray>
# include <vector>


struct glp_prob;


namespace fm
{
    class System;
    class Vector;
    typedef std::vector<Vector> Matrix;

    typedef long Value;
    typedef std::valarray<Value> ValArray;

    template <class T>
        using P = std::shared_ptr<T>;

    typedef P<void> ScopeGuard;


    struct SolveToCallback;
    struct EliminateCallback;
    struct MinimizeCallback;


    class matrix_parse_error : public std::runtime_error
    {
    public:
        using runtime_error::runtime_error;
    };

    class matrix_size_error : public std::runtime_error
    {
    public:
        using runtime_error::runtime_error;
    };


    // Minimization problem
    class Problem
    {
        P<glp_prob> prob;

        void set_mat_row(int i, const Vector&);
    public:
        size_t num_cols;

        Problem();
        explicit Problem(size_t num_cols);

        void add_equality(const Vector&);
        void add_inequality(const Vector&);
        void del_row(int i);

        bool is_redundant(const Vector&) const;
        bool dual(const Vector&, std::vector<double>&) const;
    };


    class System
    {
        // Unique ID for the current row (used for redundancy check short-cut):
        size_t next_row_id = 0;

        int get_rank(int) const;
    public:
        Matrix ineqs;
        size_t num_cols;

        explicit System(size_t nb_lines, size_t nb_cols);
        System(Matrix);

        System(System&&) = default;
        System& operator = (System&&) = default;

        // no copying
        System(const System&) = delete;
        System& operator = (const System&) = delete;

        // use this instead
        System copy() const;

        void clear(size_t nb_lines);

        void add_inequality(Vector&& v);
        void add_equality(Vector&& v);

        Problem problem() const;

        void solve_to(const SolveToCallback&, int to);
        void eliminate(const EliminateCallback&, int i);
        void minimize(const MinimizeCallback&);

        friend std::ostream& operator << (std::ostream&, const System&);
    };


    class Vector
    {
    public:
        ValArray values;

    private:
        Vector() = default;
    public:
        explicit Vector(size_t size);

        // enable move semantics
        Vector(Vector&&) = default;
        Vector& operator = (Vector&&) = default;

        // delete copy semantics to make sure, no excessive copying will
        // occur on re-allocation of the containing vector:
    private:
        Vector(const Vector&) = default;
        Vector& operator = (const Vector&) = delete;
    public:

        // use this instead
        Vector copy() const;

        bool empty() const;

        size_t size() const;
        void set(size_t i, Value n);
        Value get(size_t i) const;

        Vector eliminate(const Vector& v, size_t i) const;
        void remove(size_t i);
        void normalize();

        Vector injection(size_t dim, size_t shift=0) const;

        friend std::ostream& operator << (std::ostream&, const Vector&);
        friend Vector scaled_addition(const Vector& v0, Value s0,
                                      const Vector& v1, Value s1);
        friend bool operator == (const Vector&, const Vector&);
    };

    Vector scaled_addition(const Vector& v0, Value s0,
                           const Vector& v1, Value s1);
    ValArray scaled_addition(const ValArray& v0, Value s0,
                             const ValArray& v1, Value s1);

    size_t num_elemental_inequalities(size_t num_vars);
    fm::System elemental_inequalities(size_t num_vars);
    void set_initial_state_iid(fm::System& s, size_t width);
    void add_causal_constraints(fm::System& s, size_t width);

    int get_num_cols(const Matrix& matrix);
    int get_num_vars(const Matrix& matrix);
    Matrix copy_matrix(const Matrix& m);
    Problem problem(const Matrix& m, int num_vars);
    Matrix minimize_system(const Matrix& sys);

    fm::Vector parse_vector(std::string line);
    Matrix parse_matrix(const std::vector<std::string>& lines);


    // status/control callbacks

    struct CallbackBase {
        virtual ~CallbackBase() {}
    };

    struct SolveToCallback : CallbackBase {
        virtual ScopeGuard start_step(int step) const;
        virtual P<EliminateCallback> start_eliminate(int index) const;
    };

    struct EliminateCallback : CallbackBase {
        virtual ScopeGuard start_append(int z, int p, int n) const;
        virtual ScopeGuard start_check(int index) const;
    };

    struct MinimizeCallback : CallbackBase {
        virtual ScopeGuard start_round(int i) const;
    };

    struct StatusOutput
    {
        std::ostream* o;
        System* s;
    };

    struct SolveToStatusOutput : SolveToCallback, StatusOutput
    {
        SolveToStatusOutput(std::ostream&, System&, int to);
        ~SolveToStatusOutput();
        P<EliminateCallback> start_eliminate(int index) const;
    };

    struct EliminateStatusOutput : EliminateCallback, StatusOutput
    {
        EliminateStatusOutput(std::ostream&, System&);
        ~EliminateStatusOutput();
        ScopeGuard start_append(int z, int p, int n) const;
    };

    struct MinimizeStatusOutput : MinimizeCallback, StatusOutput
    {
        int num_orig;

        MinimizeStatusOutput(std::ostream&, System&);
        ~MinimizeStatusOutput();
        ScopeGuard start_round(int i) const;
    };

}

#endif  // include guard
