// Enumerate information inequalities using a Monte-carlo method:
//
// First, randomly drop a specified number of inequalities from the original
// set, then eliminate using the reduced set.
//
// Repeat this multiple times until the full system (=reference) is
// discovered.

#include <cmath>            // sqrt
#include <cstdlib>          // atol
#include <cstddef>
#include <fstream>

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


using boost::chrono::nanoseconds;
typedef boost::chrono::duration<double> seconds;
using boost::timer::cpu_timer;


struct Timeout
{
    cpu_timer timer;
    seconds limit;

    Timeout(seconds l)
        : limit(l)
    {
    }

    seconds elapsed() const
    {
        auto elapsed = timer.elapsed();
        return nanoseconds(elapsed.user + elapsed.system);
    }

    bool operator() () const
    {
        return elapsed() > limit;
    }
};


struct SolveToTimelimit : fm::SolveToStatusOutput
{
    Timeout timeout;

    typedef fm::SolveToStatusOutput super;

    SolveToTimelimit(const fm::IO& io, seconds timelimit)
        : super(io)
        , timeout(timelimit)
    {
    }

    fm::SG start_step(int step) const                   override
    {
        if (timeout()) {
            throw timeout_error("Running out of time.");
        }
        return super::start_step(step);
    }
    
};


fm::Matrix random_elimination(fm::System system, int num_drop,
                              seconds timelimit, const fm::IO& io)
{
    fm::Matrix& matrix = system.ineqs;
    int num_vars = matrix[0].size();
    int width = intlog2(num_vars)/2;
    int solve_to = 1<<width;

    std::random_device rd;
    std::default_random_engine random_engine(rd());
    while (num_drop-- && !matrix.empty()) {
        std::uniform_int_distribution<int> random_dist(0, matrix.size());
        int index = random_dist(random_engine);
        matrix.erase(matrix.begin() + index);
    }

    fm::solve_to{system, solve_to}.run(SolveToTimelimit(io, timelimit));
    fm::minimize{system}.run(fm::MinimizeStatusOutput(io));
    return move(system.ineqs);
}


void merge(fm::System& s, fm::Matrix& m, const fm::IO& io)
{
    fm::Problem lp = s.problem();
    for (auto&& v : m) {
        if (!lp.is_redundant(v)) {
            lp.add_inequality(v);
            s.add_inequality(move(v));
        }
    }
    m.clear();
    fm::minimize{s}.run(fm::MinimizeStatusOutput(io));
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


int main(int argc, char** argv, char** env)
try
{
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " NUM_DROP INITIAL ACCUM" << endl;
        return 1;
    }
    util::AutogenNotice gen(argc, argv);

    fm::IO io(&cerr);

    int num_drop = atol(argv[1]);;
    fm::System init_state = fm::parse_matrix(util::read_file(argv[2]));
    seconds timelimit(5*60);

    fm::Matrix result = random_elimination(move(init_state), num_drop,
                                           timelimit, io);
    fm::System accum = fm::parse_matrix(util::read_file(argv[3]));
    merge(accum, result, io);

    ofstream out(argv[3]);
    out << gen.str() << endl;
    out << accum << endl;

    return 0;
}
catch (...)
{
    throw;
}
