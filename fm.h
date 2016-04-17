// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#ifndef __FM_H__INCLUDED__
#define __FM_H__INCLUDED__

# include <iostream>
# include <memory>      // shared_ptr
# include <valarray>
# include <vector>

# include "lp.h"
# include "linalg.h"

# define EMPTY(type) { return type(); }


// External
namespace terminal { class Input; }


// Local

namespace fm
{
    // import
    template <class T> using P = std::shared_ptr<T>;
    template <class T> using Vec = la::Vector<T>;
    template <class T> using Mat = la::Matrix<T>;
    using lp::Problem;

    // helper
    typedef P<void> ScopeGuard, SG;

    // export
    class System;
    class Vector;
    typedef std::vector<Vector> Matrix;

    typedef int Value;
    typedef Vec<Value> ValArray;


    struct SolveToCallback;
    struct EliminateCallback;
    struct MinimizeCallback;


    class System
    {
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
        Vector(ValArray);

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

    size_t num_elemental_inequalities(size_t num_vars);
    fm::System elemental_inequalities(size_t num_vars);
    void set_initial_state_iid(fm::System& s, size_t width, size_t offset);
    void add_causal_constraints(fm::System& s, size_t width, size_t branches=2, bool cyclic=true);

    int get_num_cols(const Matrix& matrix);
    int get_num_vars(const Matrix& matrix);
    Matrix copy_matrix(const Matrix& m);
    Problem problem(const Matrix& m, int num_vars);
    Matrix minimize_system(const Matrix& sys);


    Vector parse_vector(std::string line);
    Matrix parse_matrix(const std::vector<std::string>& lines);

    // status/control callbacks

    struct CallbackBase {
        virtual ~CallbackBase() {}
    };

    struct minimize
    {
        System& sys;

        struct Callback : CallbackBase {
            virtual SG enter(minimize*) const EMPTY(SG);
            virtual SG start_round(int i) const EMPTY(SG);
        };
        void run(const Callback& cb=Callback());
    };

    struct eliminate
    {
        System& sys;
        int index;

        struct Callback : CallbackBase {
            virtual SG enter(eliminate*) const EMPTY(SG);
            virtual SG start_append(int z, int p, int n) const EMPTY(SG);
            virtual SG start_check(int index) const EMPTY(SG);
        };
        void run(const Callback& cb=Callback());
    };

    typedef P<eliminate::Callback> EliminatePtr;

    struct solve_to
    {
        System& sys;
        int to;
        int get_rank(int) const;

        struct Callback : CallbackBase {
            virtual SG enter(solve_to*) const EMPTY(SG);
            virtual SG start_step(int step) const EMPTY(SG);
            virtual EliminatePtr start_eliminate(int index) const;
        };
        void run(const Callback& cb=Callback());
    };

    typedef P<terminal::Input> InputPtr;

    struct IO
    {
        std::ostream* out;
        InputPtr inp;
        IO(std::ostream*, InputPtr=InputPtr());
    };

    struct MinimizeStatusOutput : minimize::Callback, IO
    {
        mutable System* sys;
        mutable int num_orig;
        MinimizeStatusOutput(IO io) : IO(io) {}
        ~MinimizeStatusOutput();
        SG enter(minimize*) const                       override;
        SG start_round(int i) const                     override;
    };

    struct EliminateStatusOutput : eliminate::Callback, IO
    {
        mutable System* sys;
        EliminateStatusOutput(IO io) : IO(io) {}
        ~EliminateStatusOutput();
        SG enter(eliminate*) const                      override;
        SG start_append(int z, int p, int n) const      override;
    };

    struct SolveToStatusOutput : solve_to::Callback, IO
    {
        mutable System* sys;
        SolveToStatusOutput(IO io) : IO(io) {}
        ~SolveToStatusOutput();
        SG enter(solve_to*) const                       override;
        SG start_step(int step) const                   override;
        EliminatePtr start_eliminate(int index) const   override;
    };

}

#endif  // include guard
