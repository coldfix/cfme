// Check if two systems of inequalities are equivalent
#include "fm.h"

#include <fstream>
#include <string>
#include <vector>

#include "util.h"
#include "error.h"

using namespace std;
using fm::Matrix;


//----------------------------------------
// business logic
//----------------------------------------

Matrix unimplied(const Matrix& a, const Matrix& b)
{
    Matrix r;
    if (a.empty() && b.empty())
        return r;
    int num_vars_a = get_num_vars(a);
    int num_vars_b = get_num_vars(b);
    int num_vars;
    if (num_vars_a == -1) {
        num_vars = num_vars_b;
    }
    else if (num_vars_b == -1) {
        num_vars = num_vars_a;
    }
    else {
        _assert<la::size_error>(num_vars_a == num_vars_b,
                "systems must not differ in size");
        num_vars = num_vars_a;
    }
    fm::Problem lp = problem(a, num_vars);
    for (auto&& v : b) {
        if (!lp.is_redundant(v)) {
            r.push_back(v.copy());
        }
    }
    return r;
}


bool check_implies(string label_a, const Matrix& sys_a,
                   string label_b, const Matrix& sys_b)
{
    Matrix missing = unimplied(sys_a, sys_b);
    if (missing.empty()) {
        cout << label_a << " implies " << label_b << endl;
        return true;
    }
    cout << label_a << " misses the following parts of " << label_b << ":"
        << endl;
    for (auto&& v : missing) {
        cout << "  " << v << endl;
    }
    return false;
}


int main(int argc, char** argv, char** env)
try
{
    int error_level = 0;

    if (argc == 3) {
        string file_a = argv[1];
        string file_b = argv[2];
        Matrix sys_a = fm::parse_matrix(util::read_file(file_a));
        Matrix sys_b = fm::parse_matrix(util::read_file(file_b));
        string label_a = "A";
        string label_b = "B";
        if (!check_implies(label_a, sys_a, label_b, sys_b)) {
            error_level |= 1;
        }
        if (!check_implies(label_b, sys_b, label_a, sys_a)) {
            error_level |= 2;
        }
    }

    else {
        cout << "Usage: check FILENAME FILENAME" << endl;
        error_level = 4;
    }

    return error_level;
}
catch (...)
{
    throw;
}
