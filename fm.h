// Basic Fourier-Motzkin C++ API (eliminates variables from a system of
// inequalities).

#ifndef __FM_H__INCLUDED__
#define __FM_H__INCLUDED__

# include <memory>          // shared_ptr


// forward declarations
struct s_fm_system;
struct s_fm_vector;


namespace fm
{
    template <class T>
        using P = std::shared_ptr<T>;


    class System;
    class Vector;


    class System
    {
        P<s_fm_system> sys;

    public:
        explicit System(s_fm_system* system);

        static System create(size_t nb_lines, size_t nb_cols);

        Vector row(int row);
        System solution_to(int to);
        void print(FILE* stream=NULL);
    };

    class Vector
    {
        s_fm_vector* vec;

    public:
        explicit Vector(s_fm_vector* vector);
        ~Vector();

        void set(int i, int n);
    };
}

#endif  // include guard
