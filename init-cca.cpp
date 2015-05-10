// - initialize a system of 2*W random variables with elemental inequalities
// - add the structural constraints of a two layered periodic CCA of width=W
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
    fm::add_causal_constraints(system, width);
    fm::Matrix& ineqs = system.ineqs;
    size_t num_orig = ineqs.size();

    fm::Problem lp = system.problem();
    for (int i = ineqs.size()-1; i >= 0; --i) {
        lp.del_row(i+1);
        if (lp.is_redundant(ineqs[i])) {
            ineqs.erase(ineqs.begin() + i);
        }
        else {
            lp.add_inequality(ineqs[i]);
        }

        terminal::clear_current_line(cerr);
        cerr << "Minimizing: " << num_orig << " -> " << ineqs.size()
            << "  (i=" << i << ")"
            << flush;
    }
    terminal::clear_current_line(cerr);
    cerr << "Minimizing: " << num_orig << " -> " << ineqs.size()
        << " (DONE)"
        << endl;

    cout << gen.str() << endl;
    cout << system << endl;
}
catch (...)
{
    throw;
}
