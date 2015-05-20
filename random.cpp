// Enumerate information inequalities using a Monte-carlo method:
//
// First, randomly drop a specified number of inequalities from the original
// set, then eliminate using the reduced set.
//
// Repeat this multiple times until the full system (=reference) is
// discovered.

#include <cstdlib>          // atol
#include <cstddef>
#include <iomanip>          // setw
#include <chrono>

#include <random>

#include <boost/timer/timer.hpp>

#include "fm.h"
#include "util.h"
#include "number.h"

using namespace std;


struct timeout_error : public std::runtime_error
{
    using runtime_error::runtime_error;
};


typedef chrono::steady_clock Clock;
typedef chrono::time_point<Clock> Time;
using std::chrono::seconds;


struct Timeout
{
    Time start;
    seconds limit;

    Timeout(seconds l)
        : start(Clock::now())
        , limit(l)
    {
    }

    bool operator() () const
    {
        return Clock::now() - start > limit;
    }
};


struct SolveToTimelimit : fm::solve_to::Callback
{
    Timeout timeout;

    SolveToTimelimit(seconds timelimit)
        : timeout(timelimit)
    {
    }

    fm::SG start_step(int step) const                   override
    {
        if (timeout()) {
            throw timeout_error("Running out of time.");
        }
        return fm::SG();
    }
};


fm::Matrix random_elimination(fm::System system, int num_drop)
{
    fm::Matrix& matrix = system.ineqs;
    size_t num_vars = matrix[0].size();
    size_t width = intlog2(num_vars)/2;
    size_t solve_to = 1<<width;

    std::random_device rd;
    std::default_random_engine random_engine(rd());
    while (num_drop-- && !matrix.empty()) {
        std::uniform_int_distribution<int> random_dist(0, matrix.size());
        int index = random_dist(random_engine);
        matrix.erase(matrix.begin() + index);
    }

    fm::solve_to{system, solve_to}.run(SolveToTimelimit(seconds(30)));
    fm::minimize{system}.run();
    return move(system.ineqs);
}


int count_nontrivial(const fm::System& a, const fm::Matrix& b)
{
    int count = 0;
    fm::Problem lp = a.problem();
    for (auto&& v : b) {
        count += !lp.is_redundant(v);
    }
    return count;
}


struct Result
{
    int number_of_steps = 0;
    std::vector<int> num_found;
    std::vector<int> num_nontrivial;
    std::vector<int> missing_nontrivial;
    Timeout timeout;

    int width;
    fm::System discovery;
    fm::System ref_solution;
    fm::System elemental;
    bool finished = false;
    int num_timeouts = 0;
    int num_nontriv;
    int num_missing;

    Result(fm::System ref)
        : width(intlog2(ref.num_cols))
        , discovery(ref.ineqs.size(), ref.num_cols)
        , ref_solution(move(ref))
        , elemental(fm::elemental_inequalities(width))
        , timeout(seconds(60))
    {
    }

    void add(fm::Matrix m)
    {
        for (auto&& v : m) {
            discovery.add_inequality(v.copy());
        }
        fm::minimize{discovery}.run();

        num_nontriv = count_nontrivial(elemental, m);
        num_missing = count_nontrivial(discovery, ref_solution.ineqs);
        finished = num_missing == 0;

        ++number_of_steps;
        num_found.push_back(m.size());
        num_nontrivial.push_back(num_nontriv);
        missing_nontrivial.push_back(num_missing);
        cerr << "i=" << number_of_steps-1
            << ", found " << setw(2) << m.size()
            << ", total " << setw(2) << discovery.ineqs.size()
            << endl;
    }

    void run(fm::System init_state, int num_drop)
    {
        while (!finished && !timeout()) {
            try {
                add(random_elimination(init_state.copy(), num_drop));
            }
            catch (timeout_error& e) {
                ++num_timeouts;
            }
        }
    }
};



int main(int argc, char** argv, char** env)
try
{
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " NUM_DROP INITIAL REFERENCE" << endl;
        return 1;
    }
    util::AutogenNotice gen(argc, argv);

    size_t num_drop = atol(argv[1]);;
    fm::System init_state = fm::parse_matrix(util::read_file(argv[2]));
    fm::System ref_solution = fm::parse_matrix(util::read_file(argv[3]));

    size_t width = intlog2(init_state.num_cols)/2;
    size_t solve_to = 1<<width;

    Result r(ref_solution.copy());

    r.run(init_state.copy(), num_drop);

    cout << gen.str() << endl;
    return 0;
}
catch (...)
{
    throw;
}
