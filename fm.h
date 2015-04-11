// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#ifndef __FM_H__INCLUDED__
#define __FM_H__INCLUDED__

# include <iostream>
# include <set>
# include <valarray>
# include <vector>


namespace fm
{
    class System;
    class Vector;

    typedef long Value;


    class System
    {
    public:
        std::vector<Vector> ineqs;
        std::vector<Vector> eqns;
        size_t num_cols;

        explicit System(size_t nb_lines, size_t nb_cols);

        void clear();

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
