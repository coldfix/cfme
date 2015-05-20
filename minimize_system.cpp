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



int main(int argc, char** argv, char** env)
try
{
    using namespace std;

    fm::System system = fm::parse_matrix(util::read_file(cin));

    util::AutogenNotice gen(argc, argv);
    fm::minimize{system}.run(fm::MinimizeStatusOutput(&cerr));

    cout << gen.str() << endl;
    cout << system << endl;
}
catch (...)
{
    throw;
}
