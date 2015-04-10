// Thin C++ wrapper classes around the API of the FM library.

#include <fm/system.h>
#include <fm/solution.h>
#include <fm/vector.h>

#include <fm/solver.h>


#include "fm.h"


namespace fm
{

    // class System

    System::System(s_fm_system_t* system)
        : sys(system, fm_system_free)
    {
    }

    System System::create(size_t nb_lines, size_t nb_cols)
    {
        return System(fm_system_alloc(nb_lines, nb_cols));
    }

    Vector System::row(int row)
    {
        return Vector(sys->lines[row]);
    }

    System System::solution_to(int to)
    {
        // Solve the given system. Use FM_SOLVER_AUTO_SIMPLIFY to help the
        // solver to scale automatically
        int flags=FM_SOLVER_FAST;
        P<s_fm_solution> solution(
                fm_solver_solution_to(sys.get(), flags, to),
                fm_solution_free);
        return System(fm_solution_to_system(solution.get()));
    }

    void System::print(FILE* stream)
    {
        if (!stream)
            stream = stdout;
        fm_system_print(stream, sys.get());
    }


    // class Vector

    Vector::Vector(s_fm_vector_t* vector)
        : vec(vector)
    {
        // this is an inequality
        set(0, 1);
    }

    Vector::~Vector()
    {
        fm_vector_compute_key(&vec->key, vec);
    }

    void Vector::set(int i, int n)
    {
        vec->vector[i].num = n;
    }

}

// vim: path=.,/usr/include,/opt/fm-lib/include
