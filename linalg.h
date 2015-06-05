#ifndef __LINALG_H__INCLUDED__
#define __LINALG_H__INCLUDED__

# include <algorithm>   // copy
# include <iterator>    // istream_iterator, back_inserter
# include <iostream>
# include <string>
# include <sstream>     // istringstream
# include <valarray>
# include <vector>

# include "error.h"     // _assert, parse_error, size_error


namespace la
{

    template <class T> using Vector = std::valarray<T>;
    template <class T> using Matrix = std::vector<Vector<T>>;

    template <class T>
    Vector<T> parse_vector(std::string line)
    {
        typedef std::istream_iterator<T> iit;
        assert_eq(line.front(), '[', parse_error, "expecting '['", line);
        assert_eq(line.back(), ']', parse_error, "expecting ']'", line);
        line = util::trim(line.substr(1, line.size()-2));
        std::istringstream in(line);
        std::vector<T> vals;
        std::copy(iit(in), iit(), std::back_inserter(vals));
        Vector<T> r(vals.size());
        std::copy(vals.begin(), vals.end(), std::begin(r));
        return r;
    }

    template <class T>
    Matrix<T> parse_matrix(const std::vector<std::string>& lines)
    {
        Matrix<T> r;
        for (std::string line : lines) {
            line = util::remove_comment(line);
            line = util::trim(line);
            if (line.empty())
                continue;
            r.push_back(parse_vector<T>(line));
        }
        return r;
    }

    template <class S, class T>
    Vector<S> convert(const Vector<T>& v)
    {
        Vector<S> r(v.size());
        for (int i = 0; i < v.size(); ++i) {
            r[i] = v[i];
        }
        return r;
    }

    template <class T>
    Vector<T> basis_vector(int size, int axis, T value=1)
    {
        Vector<T> v(size);
        v[axis] = value;
        return v;
    }

    template <class T>
    std::ostream& print_vector(std::ostream& out, const Vector<T>& vec)
    {
        out << "[ ";
        for (auto val : vec) {
            out.width(3);
            out << val << ' ';
        }
        out << "]";
        return out;
    }

    template <class T>
    std::ostream& print_matrix(std::ostream& o, const Matrix<T>& M)
    {
        for (auto&& v : M) {
            o << v << '\n';
        }
        return o;
    }

    template <class T>
    Vector<T> scaled_addition(const Vector<T>& v0, T s0,
                              const Vector<T>& v1, T s1)
    {
        assert_eq_size(v0.size(), v1.size());
        return v0 * s0 + v1 * s1;
    }

    template <class T>
    bool equal(const Vector<T>& a, const Vector<T>& b)
    {
        assert_eq_size(a.size(), b.size());
        for (int i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }

    template <class T>
    int num_rows(const Matrix<T>& M)
    {
        return M.size();
    }

    template <class T>
    int num_cols(const Matrix<T>& M)
    {
        _assert(!M.empty(), size_error);
        return M[0].size();
    }

    template <class T>
    Matrix<T> transpose(Matrix<T> M)
    {
        int nr = num_rows(M);
        int nc = num_cols(M);
        Matrix<T> r(nc, Vector<T>(nr));
        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nc; ++j) {
                r[j][i] = M[i][j];
            }
        }
        return r;
    }

    // r = Mv
    template <class T>
    Vector<T> multiply(const Matrix<T>& M, const Vector<T>& v)
    {
        int nr = num_rows(M);
        int nc = num_cols(M);
        assert_eq_size(nc, v.size());
        Vector<T> r(nr);
        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nc; ++j) {
                r[i] += M[i][j] * v[j];
            }
        }
        return r;
    }

    // r = vM
    template <class T>
    Vector<T> multiply(const Vector<T>& v, const Matrix<T>& M)
    {
        int nr = num_rows(M);
        int nc = num_cols(M);
        assert_eq_size(nr, v.size());
        Vector<T> r(nc);
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nr; ++j) {
                r[i] +=  v[j] * M[j][i];
            }
        }
        return r;
    }

    template <class T>
    Vector<T> embed(const Vector<T>& v, int dim, int shift=0)
    {
        _assert(dim >= v.size()+shift, size_error);
        Vector<T> r(dim);
        for (int i = 0; i < v.size(); ++i) {
            r[i+shift] = v[i];
        }
        return r;
    }


    template <class T>
    Vector<T> project(const Vector<T>& v, int dim_to, int dim_from=0)
    {
        Vector<T> r(dim_to - dim_from);
        for (int i = 0; i < r.size(); ++i) {
            r[i] = v[i+dim_from];
        }
        return r;
    }

}

#endif // include guard
