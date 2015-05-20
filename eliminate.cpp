// - read a system of inequalities from STDIN
// - eliminate columns from the right until reaching size=SOLVE_TO
// - minimize resulting system of inequalities
// - print system to STDOUT
//
// Status updates are shown on STDERR

#include <cstdlib>          // atol
#include <cstddef>
#include <iomanip>          // setw
#include <iostream>
#include <utility>          // move

#include "fm.h"
#include "util.h"
#include "number.h"         // intlog2


using std::cout;
using std::cerr;
using std::endl;
using std::vector;


struct RecordOrder : fm::SolveToStatusOutput
{
    vector<int>* recorded_order;

    typedef fm::SolveToStatusOutput super;

    RecordOrder(const fm::IO& io, vector<int>* r)
        : super(io)
        , recorded_order(r)
    {
    }

    fm::EliminatePtr start_eliminate(int index) const override
    {
        recorded_order->push_back(index);
        return super::start_eliminate(index);
    }
};


int main(int argc, char** argv, char** env)
try
{
    util::AutogenNotice gen(argc, argv);

    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " SOLVE_TO" << endl;
        return 1;
    }

    size_t solve_to = std::atol(argv[1]);
    size_t width = intlog2(solve_to);

    fm::IO io(&cerr);

    fm::System system = fm::parse_matrix(util::read_file(std::cin));

    // make a copy that can be used later to verify that inequalities
    // are indeed implied (consistency check for FM algorithm):
    fm::Problem orig_lp = system.problem();

    vector<int> recorded_order;
    fm::solve_to{system, solve_to}.run(RecordOrder(io, &recorded_order));
    fm::minimize{system}.run(fm::MinimizeStatusOutput(io));

    cerr << "Reduced to "
        << system.ineqs.size() << " inequalities and "
        << "Expecting " << fm::num_elemental_inequalities(width)
        << " elemental inequalities.\n"
        << endl;

    // used to remove inequalities implied by elemental inequalities on
    // the reduced space:
    fm::System target = fm::elemental_inequalities(width);

    // consistency checks
    cerr << "Perform consistency checks: " << endl;
    cerr << " - Search for false positives" << endl;
    bool consistent = true;
    for (auto&& v : system.ineqs) {
        if (!orig_lp.is_redundant(v.injection(orig_lp.num_cols))) {
            cerr << "   FALSE: " << v << endl;
            consistent = false;
        }
    }
    cerr << " - Search for undiscovered elemental inequalities" << endl;
    fm::Problem sys_prob = system.problem();
    for (auto&& v : target.ineqs) {
        if (!sys_prob.is_redundant(v)) {
            cerr << "   UNDISCOVERED: " << v << endl;
            consistent = false;
        }
    }
    cerr << endl;
    if (!consistent) {
        return 1;
    }

    // enumerate non-trivial constraints
    fm::Matrix non_trivial;
    cerr << "Filtering non-trivial inequalities." << endl;
    fm::Problem tgt_prob = target.problem();
    for (auto&& v : system.ineqs) {
        if (tgt_prob.is_redundant(v))
            continue;
        tgt_prob.add_inequality(v.copy());
        non_trivial.push_back(v.copy());
    }
    cerr << endl;

    cout << gen.str() << endl;
    cout << "\n# Elimination order:";
    for (int i = 0; i < recorded_order.size(); ++i) {
        if (i % 10 == 0) {
            cout << "\n#   ";
        }
        cout << ' ' << std::setw(3) << recorded_order[i];
    }
    cout << "\n" << endl;

    system.ineqs = move(non_trivial);
    cout << system << endl;

    return 0;
}
catch (...)
{
    throw;
}

