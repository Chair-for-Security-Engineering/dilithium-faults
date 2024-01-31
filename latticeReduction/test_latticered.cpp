#include "lattice_red.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <exception>
#include <cstring>
#include <cstdlib>
#include <fplll.h>


using namespace std;
using namespace NTL;
using namespace fplll;


int main()
{
    init_ntl();

    int nr_known = 240;
    int block_size = 15;
    int success = 0;
    cout << "Number of known NTT coefficients per s1 component " << nr_known << endl;
    cout << "Block size " << block_size << endl;
    for (auto i = 0; i < 2; i++)
    {
        success |= test_latticered(block_size, nr_known);
    }

    nr_known = 128;
    block_size = 30;
    success = 0;
    cout << "Number of known NTT coefficients per s1 component " << nr_known << endl;
    cout << "Block size " << block_size << endl;
    for (auto i = 0; i < 2; i++)
    {
        success |= test_latticered(block_size, nr_known);
    }

    return 0;
}