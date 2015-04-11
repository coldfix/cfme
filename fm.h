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
    public:
        std::vector<Vector> ineqs;
        std::vector<Vector> eqns;
        size_t num_cols;

        explicit System(size_t nb_lines, size_t nb_cols);

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
        std::valarray<Value> values;
        std::set<size_t> comb;

    private:
        Vector() = default;
    public:
        explicit Vector(size_t size, int id=-1);

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

        friend std::ostream& operator << (std::ostream&, const Vector&);
    };

}

#endif  // include guard
