// - initialize a system of 2*W random variables with elemental inequalities
// - treat this system as a two layered periodic CCA of width W and add the
//   corresponding causal constraints
// - set the first layer to be mutual independent
// - minimize the system of inequalites
// - print all vectors to STDOUT

#include <cstdlib>      // atol
#include <iostream>
#include "fm.h"

#include "util.h"

using namespace std;


int usage(int argc, char** argv)
{
    cerr << "Usage: " << argv[0] << " WIDTH [NUM_LINKS [--flat | --cyclic]" << endl;
    return 1;
}


int main(int argc, char** argv, char** env)
try
{
    if (argc < 2 || argc > 4) {
        return usage(argc, argv);
    }

    size_t width = atol(argv[1]);
    size_t links = argc >= 3 ?  atol(argv[2]) : 2;
    bool cyclic = true;

    if (argc >= 4) {
        if (argv[3] == string("--flat")) {
            cyclic = false;
        }
        else if (argv[3] == string("--cyclic")) {
            cyclic = true;
        }
        else {
            return usage(argc, argv);
        }
    }

    size_t num_final = cyclic ? width : width-1;
    size_t num_vars = width + num_final;

    util::AutogenNotice gen(argc, argv);

    fm::System system = fm::elemental_inequalities(num_vars);

    fm::set_initial_state_iid(system, width, num_final);
    fm::add_causal_constraints(system, width, links, cyclic);
    fm::minimize{system}.run(fm::MinimizeStatusOutput(&cerr));

    cout << gen.str() << endl;
    cout << system << endl;
}
catch (...)
{
    throw;
}
