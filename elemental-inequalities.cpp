#include <iostream>

#include "fm.h"

using namespace std;


int main(int argc, char** argv)
{
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " NUM_VARS" << endl;
        return 1;
    }

    size_t num_vars = atol(argv[1]);
    fm::System eli = fm::elemental_inequalities(num_vars);
    cout << eli << endl;
    return 0;
}
