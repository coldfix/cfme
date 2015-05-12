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


fm::Vector shifted(const fm::Vector& vec, int width, int shift)
{
    fm::Vector res(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        res.set(shifted(i, width, shift), vec.get(i));
    }
    return res;
}


bool check_shift_invariance(const Matrix& sys)
{
    int num_vars = get_num_vars(sys);
    int width = num_vars/2;

    vector<bool> checked(sys.size());

    bool success = true;
    for (int i = 0; i < sys.size(); ++i) {
        if (checked[i]) {
            continue;
        }

        std::vector<int> missing_shifts;

        for (int shift = 1; shift < width; ++shift) {
            fm::Vector vec = shifted(sys[i], width, shift);

            if (vec == sys[i]) {
                continue;
            }

            bool found = false;
            for (int j = i+1; j < sys.size(); ++j) {
                if (sys[j] == vec) {
                    found = true;
                    checked[j] = true;
                    break;
                }
            }
            if (!found) {
                missing_shifts.push_back(shift);
                success = false;
            }
            if (!missing_shifts.empty()) {
                cerr << "For vector: " << sys[i] << "\n";

                for (auto shift : missing_shifts) {
                    cerr << "  no shift: " 
                        << shifted(sys[i], width, shift)
                        << " (shift=" << shift << ")"
                        << endl;
                }
            }
        }

        checked[i] = true;
    }
    return success;
}


int main(int argc, char** argv, char** env)
try
{
    int error_level = 0;

    if (argc == 2) {
        string file = argv[1];
        Matrix sys = fm::parse_matrix(util::read_file(file));
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
