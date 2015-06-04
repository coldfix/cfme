// Check if two systems of inequalities are equivalent
#include "fm.h"

#include <fstream>
#include <string>
#include <vector>

#include "util.h"
#include "error.h"
#include "number.h"

using namespace std;
using fm::Matrix;
using fm::System;


fm::Vector shifted(const fm::Vector& vec, int width, int shift)
{
    fm::Vector res(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        res.set(shifted(i, width, shift), vec.get(i));
    }
    return res;
}


bool check_shift_invariance(const fm::System& sys)
{
    const fm::Matrix& mat = sys.ineqs;
    int num_vars = get_num_vars(mat);
    int width = num_vars/2;

    fm::Problem lp = sys.problem();

    bool success = true;
    for (int i = 0; i < mat.size(); ++i) {
        std::vector<int> missing_shifts;

        for (int shift = 1; shift < width; ++shift) {
            fm::Vector vec = shifted(mat[i], width, shift);

            if (vec == mat[i]) {
                continue;
            }

            // The minimized system of inequalities turns out not to be
            // unique. Any of the shifted variants of an inequality might be
            // removed while still being valid. Therefore, it's insufficient
            // to just search for the vector. An LP must be solved instead:
            if (!lp.is_redundant(vec.values)) {
                missing_shifts.push_back(shift);
                success = false;
            }
            if (!missing_shifts.empty()) {
                cerr << "For vector: " << mat[i] << "\n";

                for (auto shift : missing_shifts) {
                    cerr << "  no shift: " 
                        << shifted(mat[i], width, shift)
                        << " (shift=" << shift << ")"
                        << endl;
                }
            }
        }
    }
    return success;
}


int main(int argc, char** argv, char** env)
try
{
    int error_level = 0;

    if (argc == 2) {
        string file = argv[1];
        System sys = fm::parse_matrix(util::read_file(file));
        if (!check_shift_invariance(sys)) {
            error_level = 1;
        }
    }

    else {
        cout << "Usage: check FILENAME [FILENAME]" << endl;
        error_level = 2;
    }

    return error_level;
}
catch (...)
{
    throw;
}
