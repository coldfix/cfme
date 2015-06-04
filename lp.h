#ifndef __LP_H__INCLUDED__
#define __LP_H__INCLUDED__

# include <memory>
# include "linalg.h"


struct glp_prob;


namespace lp
{
    template <class T> using P = std::shared_ptr<T>;

    typedef la::Vector<double> Vector;
    typedef la::Vector<int>   iVector;

    enum Status {
        UNDEF=1,/* solution is undefined */
        FEAS,   /* solution is feasible */
        INFEAS, /* solution is infeasible */
        NOFEAS, /* no feasible solution exists */
        OPT,    /* solution is optimal */
        UNBND,  /* solution is unbounded */
    };

    // Linear minimization problem
    class Problem
    {
        P<glp_prob> prob;

        void set_mat_row(int i, const Vector&);
    public:
        size_t num_cols;

        Problem();
        explicit Problem(size_t num_cols);

        void add_equality(const Vector&, double rhs=0);
        void add_inequality(const Vector&, double lb=0, double ub=INFINITY);
        void del_row(int i);

        bool is_redundant(const Vector&) const;
        Status simplex(const Vector&, Vector* o=nullptr) const;
        bool dual(const Vector&, Vector&) const;

        void add_equality(const iVector&, double rhs=0);
        void add_inequality(const iVector&, double lb=0, double ub=INFINITY);
        bool is_redundant(const iVector&) const;
    };

}

#endif // include guard
