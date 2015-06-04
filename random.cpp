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
#include <iomanip>          // setw
#include <sstream>

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

    Result(fm::System ref, seconds timelimit)
        : width(intlog2(ref.num_cols))
        , discovery(ref.ineqs.size(), ref.num_cols)
        , ref_solution(move(ref))
        , elemental(fm::elemental_inequalities(width))
        , timeout(timelimit)
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

        num_found.push_back(m.size());
        num_nontrivial.push_back(num_nontriv);
        missing_nontrivial.push_back(num_missing);
        terminal::clear_current_line(cerr);
        cerr << "i=" << setw(3) << number_of_steps++
            << ", found " << setw(2) << m.size()
            << ", total " << setw(2) << discovery.ineqs.size()
            << flush;
    }

    void run(fm::System init_state, int num_drop)
    {
        timeout.timer.start();
        while (!finished && !timeout()) {
            try {
                add(random_elimination(init_state.copy(), num_drop));
            }
            catch (timeout_error& e) {
                ++num_timeouts;
            }
        }
        timeout.timer.stop();
    }
};


double sq(double v)
{
    return v * v;
}


string format(double d)
{
    ostringstream out;
    out
      << std::fixed
      << std::right
      << std::setprecision(3)
      << std::setw(8)
      << d;
    return out.str();
}


struct Statistic
{
    double sum = 0;
    double sumsq = 0;
    int count = 0;

    void add(double val)
    {
        sum += val;
        sumsq += sq(val);
        ++count;
    }

    double mean() const
    {
        return sum / count;
    }

    double var() const
    {
        return sumsq - count * sq(mean());
    }

    double stddev() const
    {
        return sqrt(var() / (count-1));
    }

    double stderr() const
    {
        return stddev() / sqrt(count);
    }

    template <class V, class F>
    static Statistic accumulate(const V& vec, F func)
    {
        Statistic s;
        for (auto&& val : vec) {
            s.add(func(val));
        }
        return s;
    }

    friend ostream& operator << (ostream& out, const Statistic& s)
    {
        return out << format(s.mean())
            << " " << format(s.stderr())
            << " " << format(s.stddev());
    }
};


struct MultiRun
{
    std::vector<Result> results;

    fm::System ref_solution;
    typedef const Result& R;

    MultiRun(fm::System ref)
        : ref_solution(move(ref))
    {
    }

    void run(fm::System init_state,
            int num_drop,
            int num_turns,
            seconds timelimit)
    {
        for (int i = 0; i < num_turns; ++i) {
            cerr << i << ":" << endl;
            Result r(ref_solution.copy(), timelimit);
            r.run(init_state.copy(), num_drop);
            results.push_back(move(r));
            cerr << endl;
        }
    }

    template <class F>
    Statistic accumulate(F func) const
    {
        return Statistic::accumulate(results, func);
    }


    Statistic runtime() const
    {
        return accumulate([](R r) { return r.timeout.elapsed().count(); });
    }

    Statistic num_steps() const
    {
        return accumulate([](R r) { return r.number_of_steps; });
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

    int num_drop = atol(argv[1]);;
    fm::System init_state = fm::parse_matrix(util::read_file(argv[2]));
    fm::System ref_solution = fm::parse_matrix(util::read_file(argv[3]));
    int num_turns = 100;
    seconds timelimit(5*60);

    MultiRun r(move(ref_solution));
    r.run(init_state.copy(), num_drop, num_turns, timelimit);

    cout << gen.str() << endl;
    cout << "#"
        << setw(3) << "N"   << "    "
        << setw(8) << "t"   << " "
        << setw(8) << "err" << " "
        << setw(8) << "dev" << "    "
        << setw(8) << "n"   << " "
        << setw(8) << "err" << " "
        << setw(8) << "dev" << endl
        << setw(4) << num_drop  << "    "
        << r.runtime()          << "    "
        << r.num_steps()        << endl
        ;
    return 0;
}
catch (...)
{
    throw;
}
