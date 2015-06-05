
#include <cmath>    // NAN
#include <glpk.h>
#include "lp.h"


namespace lp
{

    // Always use zero based indices for all local variables and parameters
    // right until passing them to GLPK.

    Problem::Problem()
    {
    }

    Problem::Problem(size_t nb_cols)
        : num_cols(nb_cols)
    {
        prob.reset(glp_create_prob(), glp_delete_prob);
        glp_set_obj_dir(prob.get(), GLP_MIN);
        glp_add_cols(prob.get(), num_cols);
        for (int j = 0; j < num_cols; ++j) {
            glp_set_col_bnds(prob.get(), j+1, GLP_FR, NAN, NAN);
        }
    }

    void Problem::set_mat_row(int i, const Vector& v)
    {
        std::vector<int> indices;
        std::vector<double> values;
        indices.reserve(v.size());
        values.reserve(v.size());
        for (int j = 0; j < v.size(); ++j) {
            double val = v[j];
            if (val) {
                indices.push_back(j+1);
                values.push_back(val);
            }
        }
        glp_set_mat_row(prob.get(), i+1,
                indices.size(), indices.data()-1, values.data()-1);
    }

    void Problem::add_equality(const Vector& v, double rhs)
    {
        int i = glp_add_rows(prob.get(), 1)-1;
        glp_set_row_bnds(prob.get(), i, GLP_FX, rhs, rhs);
        set_mat_row(i, v);
    }

    void Problem::add_inequality(const Vector& v, double lb, double ub)
    {
        int i = glp_add_rows(prob.get(), 1)-1;
        if (lb == -INFINITY && ub == INFINITY) {
            glp_set_row_bnds(prob.get(), i+1, GLP_FR, NAN, NAN);
        }
        else if (lb > -INFINITY && ub == INFINITY) {
            glp_set_row_bnds(prob.get(), i+1, GLP_LO, lb, NAN);
        }
        else if (lb == -INFINITY && ub < INFINITY) {
            glp_set_row_bnds(prob.get(), i+1, GLP_LO, NAN, ub);
        }
        else {
            glp_set_row_bnds(prob.get(), i+1, GLP_DB, lb, ub);
        }
        set_mat_row(i, v);
    }

    void Problem::del_row(int i)
    {
        glp_del_rows(prob.get(), 1, (&++i)-1);
    }

    bool Problem::is_redundant(const Vector& v) const
    {
        return simplex(v) == OPT;
    }

    Status Problem::simplex(const Vector& v, Vector* o) const
    {
        assert_eq_size(v.size(), num_cols);
        for (int i = 0; i < num_cols; ++i) {
            glp_set_obj_coef(prob.get(), i+1, v[i]);
        }
        glp_std_basis(prob.get());
        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_ERR;
        int result = glp_simplex(prob.get(), &parm);
        if (result != 0) {
            throw std::runtime_error("Error in glp_simplex.");
        }

        int status = glp_get_status(prob.get());
        if (status == GLP_OPT && o) {
            for (int i = 0; i < o->size(); ++i) {
                (*o)[i] = glp_get_col_prim(prob.get(), i+1);
            }
        }
        return (Status) status;
    }

    bool Problem::dual(const Vector& v, Vector& r) const
    {
        glp_prob* lp = prob.get();
        int num_rows = glp_get_num_rows(lp);;
        assert_eq_size(v.size(), num_cols);
        assert_eq_size(r.size(), num_rows);
        for (int i = 0; i < num_cols; ++i) {
            glp_set_obj_coef(lp, i+1, v[i]);
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
        for (int i = 0; i < num_rows; ++i) {
            r[i] = glp_get_row_dual(lp, i+1);
        }
        return true;
    }

    void Problem::add_equality(const iVector& v, double rhs)
    {
        add_equality(la::convert<double>(v), rhs);
    }

    void Problem::add_inequality(const iVector& v, double lb, double ub)
    {
        add_inequality(la::convert<double>(v), lb, ub);
    }

    bool Problem::is_redundant(const iVector& v) const
    {
        return is_redundant(la::convert<double>(v));
    }

}
