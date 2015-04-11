// Thin C++ wrapper classes around the API of the FM library.


#include <cassert>
#include <utility>      // move

#include "number.h"
#include "fm.h"


using std::move;


namespace fm
{
    // class System

    System::System(size_t nb_lines, size_t nb_cols)
        : num_cols(nb_cols)
    {
        ineqs.reserve(nb_lines);
    }

    void System::clear()
    {
        ineqs.clear();
        eqns.clear();
    }

    void System::add_equality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (!vec.empty())
            ineqs.push_back(move(vec));
    }

    void System::add_inequality(Vector&& vec)
    {
        assert(vec.size() == num_cols);
        if (!vec.empty())
            eqns.push_back(move(vec));
    }

    void System::solve_to(int to)
    {
        int step_count = 0;
        for (int i = num_cols-1; i >= to; --i) {
            eliminate(i, step_count);
        }
    }

    bool System::is_redundant(const Vector& v) const
    {
        // TODO: drop if obviously redundant (multiple of other)
        // TODO: drop redundant constraints
        return false;
    }

    void System::eliminate(int index, int& step_counter)
    {
        std::cout << "elim: " << index << " " << step_counter << std::endl;

        --num_cols;
        std::vector<Vector> _ineqs = move(ineqs);
        std::vector<Vector> _eqns = move(eqns);
        clear();

        // Partition inequality constraints into (positive, negative, zero)
        // coefficient for the given index.
        std::vector<Vector> pos, neg;
        for (auto&& vec : _ineqs) {
            Value val = vec.get(index);
            if (val > 0) {
                pos.push_back(move(vec));
            }
            else if (val < 0) {
                neg.push_back(move(vec));
            }
            else if (val == 0) {
                vec.remove(index);
                add_inequality(move(vec));
            }
        }

        std::vector<Vector> eq_with;
        eq_with.reserve(_eqns.size());
        for (auto&& vec : _eqns) {
            if (vec.get(index)) {
                eq_with.push_back(move(vec));
            }
            else {
                vec.remove(index);
                add_equality(move(vec));
            }
        }

        if (!eq_with.empty()) {
            // TODO: heuristic for choosing equation?
            Vector eq = move(eq_with.back());
            eq_with.pop_back();
            for (auto&& v : eq_with) {
                add_equality(v.eliminate(eq, index));
            }
            for (auto&& v : pos) {
                add_inequality(v.eliminate(eq, index));
            }
            for (auto&& v : neg) {
                add_inequality(v.eliminate(eq, index));
            }
        }
        else {
            int step = step_counter++;

            for (auto&& p : pos) {
                for (auto&& n : neg) {
                    Vector v = p.eliminate(n, index);
                    if (v.comb.size() > 2+step)
                        continue;
                    if (!is_redundant(v))
                        add_inequality(move(v));
                }
            }
        }
    }


    // class Vector

    Vector::Vector(size_t size, int id)
        : values(size)
    {
        if (id != -1) {
            comb.insert(id);
        }
    }

    Vector Vector::copy() const
    {
        return Vector(*this);
    }

    bool Vector::empty() const
    {
        for (const Value& x : values) {
            if (x) {
                return false;
            }
        }
        return true;
    }

    size_t Vector::size() const
    {
        return values.size();
    }

    Value Vector::get(size_t i) const
    {
        return values[i];
    }

    void Vector::set(size_t i, Value v)
    {
        values[i] = v;
    }

    void Vector::remove(size_t i)
    {
        size_t s = size();
        std::valarray<Value> r(s - 1);
        for (int j = 0; j < i; ++j) {
            r[j] = values[j];
        }
        for (int j = i+1; j < s; ++j) {
            r[j-1] = values[j];
        }
        values.swap(r);
    }

    // Inplace elimination of coefficients.
    Vector Vector::eliminate(const Vector& v, size_t i) const
    {
        Vector r;

        r.comb = comb;
        for (auto x : v.comb) {
            r.comb.insert(x);
        }

        Value a = get(i);
        Value b = v.get(i);
        Value s = -sign(a*b);
        a = abs(a);
        b = abs(b);
        Value div = gcd(a, b);

        r.values = values * (b / div) + v.values * (s * (a / div));
        r.normalize();
        r.remove(i);
        return r;
    }

    // Inplace normalization of coefficients.
    void Vector::normalize()
    {
        Value div(0);
        for (const Value& x : values) {
            if (!x) {
                continue;
            }
            if (div) {
                div = gcd<Value>(div, abs(x));
            }
            else {
                div = x;
            }
            if (div == 1) {
                return;
            }
        }
        if (div > 1) {
            values /= div;
        }
    }


    // IO

    std::ostream& operator << (std::ostream& o, const System& s)
    {
        for (auto&& v : s.ineqs) {
            o << v << '\n';
        }
        return o;
    }

    std::ostream& operator << (std::ostream& o, const Vector& v)
    {
        o << "[ ";
        for (auto val : v.values) {
            o << val << ' ';
        }
        o << "]";
        return o;
    }

}
