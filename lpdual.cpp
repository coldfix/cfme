#include <iostream>
#include "linalg.h"
#include "lp.h"
#include "number.h"

using namespace std;

typedef double D;

typedef la::Vector<D> Vector;
typedef la::Matrix<D> Matrix;


int main(int argc, char** argv)
try
{
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " SYSTEM TRIAL_VECTORS" << endl;
        return 1;
    }

    Matrix matrix = la::parse_matrix<D>(util::read_file(argv[1]));
    Matrix random = la::parse_matrix<D>(util::read_file(argv[2]));

    int num_rows = la::num_rows(matrix);
    int num_cols = la::num_cols(matrix);
    int num_vars = intlog2(num_cols);
    int proj_dim = 1<<(num_vars/2);

    // setup system

    util::AutogenNotice gen(argc, argv);

    lp::Problem dlp(num_rows);
    for (int i = 0; i < num_rows; ++i) {            // 0 ≤ y ≤ 100
        dlp.add_inequality(la::basis_vector<D>(num_rows, i), 0, 100);
    }

    Matrix columns = la::transpose(matrix);
    for (int i = proj_dim; i < num_cols; ++i) {     // (y∙L)_i = 0
        dlp.add_equality(columns[i]);
    }

    for (int k = 0; k < random.size(); ++k) {
        Vector y(num_rows);
        Vector r = la::embed(random[k], num_cols);
        Vector Lr = la::multiply(matrix, r);
        // Vector Lr = random[k];
        auto status = dlp.simplex(Lr, &y);
        if (status != lp::OPT) {
            cout << "Not optimal! Status: " << status << endl;
            continue;
        }
        Vector e = la::multiply(y, matrix);
        for (D& x : e) {
            if (abs(x) < 1e-10) {
                x = 0;
            }
        }
        la::print_vector(cout, la::project(e, proj_dim)) << endl;
    }

    cout << gen.str();
    return 0;
}
catch (std::exception& e)
{
    std::cerr << "Exception: " << e.what() << std::endl;
    throw;
}
catch (...)
{
    std::cerr << "Unknown exception!" << std::endl;
    throw;
}
