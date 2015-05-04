// Enumerate information inequalities in subsequent layers of a periodic CCA
// with the given size.

#include <cstdlib>          // atol
#include <cstddef>

#include "fm.h"


// Enumerate information inequalities in second layer of a CCA of the given
// width. The initial layer is initialized to be mutually independent. The
// layout of the CCA is as described above (c.f. `add_causal_constraints`).
bool solve(size_t width)
{
    using std::cout;
    using std::endl;

    size_t num_vars = width*2;
    size_t solve_to = 1<<width;

    cout << "Initialize CCA with N=" << width << endl;
    fm::System system = fm::elemental_inequalities(num_vars);
    fm::set_initial_state_iid(system, width);
    fm::add_causal_constraints(system, width);
    cout << endl;

    for (size_t layer = 1; ; ++layer) {

        // make a copy that can be used later to verify that inequalities
        // are indeed implied (consistency check for FM algorithm):
        fm::Problem orig_lp = system.problem();

        cout << "Eliminate layer " << layer << endl;
        system.minimize();
        system.solve_to(solve_to);
        system.minimize();
        cout << endl;

        cout << "Reduced to "
            << system.ineqs.size() << " inequalities and "
            << system.eqns.size() << " equalities.\n"
            << "Expecting " << fm::num_elemental_inequalities(width)
            << " elemental inequalities.\n"
            << endl;

        // used to remove inequalities implied by elemental inequalities on
        // the reduced space:
        fm::System target = fm::elemental_inequalities(width);

        // consistency checks
        cout << "Perform consistency checks: " << endl;
        cout << " - Search for false positives" << endl;
        bool consistent = true;
        for (auto&& v : system.ineqs) {
            if (!orig_lp.is_redundant(v.injection(orig_lp.num_cols))) {
                cout << "   FALSE: " << v << endl;
                consistent = false;
            }
        }
        cout << " - Search for undiscovered elemental inequalities" << endl;
        fm::Problem sys_prob = system.problem();
        for (auto&& v : target.ineqs) {
            if (!sys_prob.is_redundant(v)) {
                cout << "   UNDISCOVERED: " << v << endl;
                consistent = false;
            }
        }
        cout << endl;
        if (!consistent) {
            return false;
        }

        // enumerate non-trivial constraints
        std::vector<fm::Vector> extra_ineqs;
        std::vector<fm::Vector> extra_eqns;
        cout << "List non-trivial inequalities: " << endl;
        fm::Problem tgt_prob = target.problem();
        for (auto&& v : system.ineqs) {
            if (tgt_prob.is_redundant(v))
                continue;
            tgt_prob.add_inequality(v.copy());
            extra_ineqs.push_back(v.copy());
            cout << v << endl;
        }
        if (extra_ineqs.empty()) {
            cout << " - None." << endl;
        }
        cout << endl;
        cout << "List equalities: " << endl;
        for (auto&& v : system.eqns) {
            extra_eqns.push_back(v.copy());
        }
        if (extra_eqns.empty()) {
            cout << " - None." << endl;
        }
        cout << endl;

        // only trivial inequalities -> return:
        if (extra_ineqs.empty() && extra_eqns.empty()) {
            break;
        }

        cout << "Initialize layer " << layer+1 << endl;
        system = fm::elemental_inequalities(num_vars);
        fm::add_causal_constraints(system, width);
        for (auto&& v : extra_ineqs) {
            system.add_inequality(v.injection(system.num_cols, width));
        }
        for (auto&& v : extra_eqns) {
            system.add_equality(v.injection(system.num_cols, width));
        }

    }

    return true;
}



int main(int argc, char** argv, char** env)
try
{
    size_t width = 2;
    if (argc >= 2)
        width = std::atol(argv[1]);
    solve(width);
    return 0;
}
catch (...)
{
    throw;
}
