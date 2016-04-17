// - initialize a system of 2*W random variables with elemental inequalities
// - treat this system as a two layered periodic CCA of width W and add the
//   corresponding causal constraints
// - read from STDIN additional constraints for the first layer
// - minimize the system of inequalites
// - print all vectors to STDOUT

#include <cstdlib>      // atol
#include <iostream>
#include "fm.h"

#include "util.h"

using namespace std;


int usage(int argc, char** argv)
{
    cerr << "Usage: " << argv[0] << " WIDTH [NUM_LINKS [NUM_INIT]]" << endl;
    return 1;
}


int main(int argc, char** argv, char** env)
try
{
    if (argc < 2 || argc > 4) {
        return usage(argc, argv);
    }

    size_t nf = atol(argv[1]);                  // num vars in final layer
    size_t nl = argc >= 3 ? atol(argv[2]) : 2;  // num links for each variable
    size_t ni = argc >= 4 ? atol(argv[3]) : nf; // num vars in initial layer
    size_t num_vars = nf + ni;

    util::AutogenNotice gen(argc, argv);

    fm::System system = fm::elemental_inequalities(num_vars);
    fm::add_causal_constraints(system, nf, ni, nl);

    for (auto&& constraint : fm::parse_matrix(util::read_file(cin))) {
        system.add_inequality(constraint.injection(system.num_cols, num_vars));
    }

    fm::minimize{system}.run(fm::MinimizeStatusOutput(&cerr));

    cout << gen.str() << endl;
    cout << system << endl;
    return 0;
}
catch (...)
{
    throw;
}
