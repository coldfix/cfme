// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#ifndef __FM_H__INCLUDED__
#define __FM_H__INCLUDED__

# include <iostream>
# include <memory>
# include <set>
# include <valarray>
# include <vector>


struct glp_prob;


namespace fm
{
    class System;
    class Vector;

    typedef long Value;
    typedef std::valarray<Value> ValArray;

    template <class T>
        using P = std::shared_ptr<T>;


    // Minimization problem
    class Problem
    {
        P<glp_prob> prob;
        size_t num_cols;

        void set_mat_row(int i, const Vector&);
    public:
        Problem();
        explicit Problem(size_t num_cols);

        void add_equality(const Vector&);
        void add_inequality(const Vector&);

        bool is_redundant(const Vector&) const;
    };


    class System
    {
        Problem problem;

        // Unique ID for the current row (used for redundancy check short-cut):
        size_t next_row_id = 0;
    public:
        std::vector<Vector> ineqs;
        std::vector<Vector> eqns;
        size_t num_cols;

        explicit System(size_t nb_lines, size_t nb_cols);

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

        void solve_to(int to);

        bool is_redundant(const Vector&) const;
        void eliminate(int i, int& s);

        friend std::ostream& operator << (std::ostream&, const System&);
    };


    class Vector
    {
    public:
        ValArray values;
        std::set<size_t> comb;

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
    };

    Vector scaled_addition(const Vector& v0, Value s0,
                           const Vector& v1, Value s1);
    ValArray scaled_addition(const ValArray& v0, Value s0,
                             const ValArray& v1, Value s1);
}

#endif  // include guard
