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


int main(int argc, char** argv, char** env)
try
{
    if (argc != 2 && argc != 3) {
        cerr << "Usage: " << argv[0] << " WIDTH [NUM_LINKS]" << endl;
        return 1;
    }
    size_t width = atol(argv[1]);
    size_t links = argc == 3 ?  atol(argv[2]) : 2;
    size_t num_vars = 2*width;

    util::AutogenNotice gen(argc, argv);

    fm::System system = fm::elemental_inequalities(num_vars);
    fm::set_initial_state_iid(system, width);
    fm::add_causal_constraints(system, width, links);
    fm::minimize{system}.run(fm::MinimizeStatusOutput(&cerr));

    cout << gen.str() << endl;
    cout << system << endl;
}
catch (...)
{
    throw;
}
