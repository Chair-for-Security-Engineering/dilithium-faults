#include "correction_attacks.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <exception>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <numeric>




#ifndef TESTS
#define TESTS 1
#endif

using namespace std;
using namespace NTL;

typedef chrono::duration<double, ratio< 1>>   sf;

int faults;

int main()
{
    /*************************************************************/
    /* this program generates statistical values:                */
    /* how long takes the attack,                                */
    /* how many faults are needed,                               */
    /* does it fail, if yes how often                            */
    /*************************************************************/

    init_ntl();
    
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    Vec<ZZ_pX> s1_recovered;

    uint8_t rho[SEEDBYTES];
    uint8_t tr[SEEDBYTES];
    uint8_t key[SEEDBYTES];
    polyveck t0;
    polyvecl s1;
    polyveck s2;
    Vec<ZZ_pX> s1_orig;

    //values for statistics
    int success = 0;
    double secs = 0.0;

    double secs_vec[TESTS] = {};
    int faults_vec[TESTS] = {};

    
    //here the first correction attack using the skipping fault
    if (SHUFFLING == 0)
    {
        for (auto i = 0; i < TESTS; i++)
        {
            //print progress 
            if (i%100 == 0)
            {
                cout << "Starting test " << i << endl;
            }

            faults = 0;

            crypto_sign_keypair(pk, sk);
            auto tic = chrono::high_resolution_clock::now();
            try
            {
                s1_recovered = recover_s1_skip(sk, pk, MLEN);
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                continue;
            }
            auto toc = chrono::high_resolution_clock::now();

            secs = sf(toc-tic).count();
            secs_vec[i] = secs;
            faults_vec[i] = faults;
            
            //get the original key
            unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
            to_ntl(s1_orig, s1);

            //compare the recovered s1 to the original s1
            if (s1_orig == s1_recovered)
            {
                success++;
            }
        }
    }
    else //shuffling on
    {
        for (auto i = 0; i < TESTS; i++)
        {
            if (i%2 == 0)
            {
                cout << "Starting test " << i << endl;
            }

            faults = 0;

            crypto_sign_keypair(pk, sk);
            auto tic = chrono::high_resolution_clock::now();
            try
            {
                s1_recovered = recover_s1_skip_shuff(sk, pk, MLEN);
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                continue;
            }
            auto toc = chrono::high_resolution_clock::now();
            
            secs = sf(toc-tic).count();
            secs_vec[i] = secs;
            faults_vec[i] = faults;
            
            //get the original key
            unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
            to_ntl(s1_orig, s1);

            //compare the recovered s1 to the original s1
            if (s1_orig == s1_recovered)
            {
                success++;
            }
        }
    }
    

    double avg_secs = 0.0;
    double max_secs = 0.0;
    double avg_faults = 0.0;
    int max_faults = 0;

    for (auto i = 0; i < TESTS; i++)
    {
        avg_secs += secs_vec[i];
        avg_faults += faults_vec[i];
        if (secs_vec[i] > max_secs)
        {
            max_secs = secs_vec[i];
        }
        if (faults_vec[i] > max_faults)
        {
            max_faults = faults_vec[i];
        }
    }
    
    avg_secs = avg_secs/((double) TESTS);
    avg_faults = ((double) avg_faults)/((double) TESTS);

    cout << "Skipping fault correction attack: Tested s1 recovery " << TESTS << " times. Succeeded " << success << " times." << endl;    
    cout << "s1 recovery took " << avg_secs << " sec on average and " << max_secs << " sec maximum." << endl;
    cout << "s1 recovery took " << avg_faults << " faults on average and " << max_faults << " faults maximum." << endl;
    
    return 0;
}