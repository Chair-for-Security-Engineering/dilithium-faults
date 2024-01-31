#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <exception>
#include <numeric>
#include <algorithm>
#include <map>
#include <chrono>
#include <random>

#include "NTL/lzz_p.h"
#include "NTL/vec_lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
#include "NTL/mat_lzz_p.h"

#include "NTL/ZZ_p.h"
#include "NTL/vec_ZZ_p.h"
#include "NTL/vec_ZZ.h"
#include "NTL/ZZ_pE.h"
#include "NTL/ZZ_pX.h"
#include "NTL/mat_ZZ_p.h"

#include "NTL/vector.h"

#include "NTL/ZZ.h"
#include "NTL/matrix.h"
#include "NTL/LLL.h"


extern "C" {
#include "dilithium/avx2/poly.h"
#include "dilithium/avx2/polyvec.h"
#include "dilithium/avx2/packing.h"
#include "dilithium/avx2/sign.h"
#include "dilithium/avx2/randombytes.h"
#include "dilithium/avx2/params.h"
#include "dilithium/avx2/fips202.h"
}

#define NSETS 10
#define NFAULTS 1
#define NFAULTS_MINIMUM 1

using namespace std;
using namespace NTL;

typedef chrono::duration<double, ratio< 1>>   sf;

void randombytes(uint8_t *out, size_t outlen)
{
    return;
}

static const size_t unpack_map[N] = {
0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57, 2, 10, 18, 26, 34, 42, 50, 58, 3, 11, 19, 27, 35, 43, 51, 59, 4, 12, 20, 28, 36, 44, 52, 60, 5, 13, 21, 29, 37, 45, 53, 61, 6, 14, 22, 30, 38, 46, 54, 62, 7, 15, 23, 31, 39, 47, 55, 63, 64, 72, 80, 88, 96, 104, 112, 120, 65, 73, 81, 89, 97, 105, 113, 121, 66, 74, 82, 90, 98, 106, 114, 122, 67, 75, 83, 91, 99, 107, 115, 123, 68, 76, 84, 92, 100, 108, 116, 124, 69, 77, 85, 93, 101, 109, 117, 125, 70, 78, 86, 94, 102, 110, 118, 126, 71, 79, 87, 95, 103, 111, 119, 127, 128, 136, 144, 152, 160, 168, 176, 184, 129, 137, 145, 153, 161, 169, 177, 185, 130, 138, 146, 154, 162, 170, 178, 186, 131, 139, 147, 155, 163, 171, 179, 187, 132, 140, 148, 156, 164, 172, 180, 188, 133, 141, 149, 157, 165, 173, 181, 189, 134, 142, 150, 158, 166, 174, 182, 190, 135, 143, 151, 159, 167, 175, 183, 191, 192, 200, 208, 216, 224, 232, 240, 248, 193, 201, 209, 217, 225, 233, 241, 249, 194, 202, 210, 218, 226, 234, 242, 250, 195, 203, 211, 219, 227, 235, 243, 251, 196, 204, 212, 220, 228, 236, 244, 252, 197, 205, 213, 221, 229, 237, 245, 253, 198, 206, 214, 222, 230, 238, 246, 254, 199, 207, 215, 223, 231, 239, 247, 255
};

static char dirname[128];



/*******************************************************************************/
/* initialize the ntl parameter for modulo reductions                          */
/*******************************************************************************/

void init_ntl() 
{
    ZZ Q_ZZ(Q);
    ZZ_p::init(Q_ZZ);
    ZZ_pX reduction_poly;
    SetCoeff(reduction_poly, 0);
    SetCoeff(reduction_poly, N);
    ZZ_pE::init(reduction_poly);
}



/*******************************************************************************/
/* convert from dilithium types to ntl types and back                          */
/*******************************************************************************/

/*poly to/from ntl_polynom (ZZ_pX) or ntl_vec (Vec<ZZ_p>)*/
void to_ntl(ZZ_pX& out, const poly& in)
{
    out.SetLength(N);
	for(auto i = 0; i < N; i++)
    {
        SetCoeff(out, i, in.coeffs[i]);
	}
}

void to_ntl(Vec<ZZ_p>& out, const poly& in)
{
	out.SetLength(N);
	for(auto i = 0; i < N; i++)
    {
        out[i] = in.coeffs[i];
	}
}

void from_ntl(poly& out, const ZZ_pX& in)
{
	for(auto i = 0; i < N; i++)
    {
        ZZ_p c;
        c = coeff(in, i);
		uint32_t tmp;
		conv(tmp, c);
		out.coeffs[i] = tmp;
	}
}

void from_ntl(poly& out, const Vec<ZZ_p>& in)
{
	for(auto i = 0; i < N; i++)
    {
        auto c = in[i];
		uint32_t tmp;
		conv(tmp, c);
		out.coeffs[i] = tmp;
	}
}


/*polyvecl to/from ntl_vec<polynom> (Vec<ZZ_pX>) of length L*/
void to_ntl(Vec<ZZ_pX>& out, const polyvecl& in)
{
    out.SetLength(L);
	for(auto i = 0; i < L; i++)
    {
		to_ntl(out[i], in.vec[i]);
	}
}

void from_ntl(polyvecl& out, const Vec<ZZ_pX>& in)
{
    for(auto i=0; i < L; i++)
    {
        from_ntl(out.vec[i], in[i]);
    }
}

/*******************************************************************************/
/* FAULT IN A (fA)                                                             */
/*                                                                             */
/* correct_sig_fA is a modified verify function                                */
/* It takes a target coefficient in A and a fA-faulted signature               */
/* as well as the delta in the NTT representation of A                         */
/* The output is the target coeff in the NTT representation of y               */
/*                                                                             */
/* recover_s1_fA(_ntt_coeff) recover s1 using the recovered NTT(y) coeff.      */
/*******************************************************************************/

int correct_sig_fA( ZZ_p &res,
                    uint8_t *sig,
                    size_t siglen,
                    const uint8_t *m,
                    size_t mlen,
                    const uint8_t *pk,
                    const polyvecl mat[K],
                    unsigned int target_row,
                    unsigned int target_col, 
                    unsigned int target_coeff, 
                    int delta)
{
    uint8_t rho[SEEDBYTES];
    uint8_t mu[CRHBYTES];
    uint8_t c[SEEDBYTES];
    poly cp;
    polyvecl z;
    polyveck t1, w1, h;
    keccak_state state;
    int foundflag = 0;
    int recovered_ntty = -1;
    int diff = 0;

    if(siglen != CRYPTO_BYTES)
    {
        #pragma omp critical
        cerr << "bad siglen\n";
        return 0;
    }
    unpack_pk(rho, &t1, pk);
    if(unpack_sig(c, &z, &h, sig))
    {
        #pragma omp critical
        cerr << "bad unpack\n";
        return 0;
    }

    /* Compute CRH(H(rho, t1), msg) */
    shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    shake256_init(&state);
    shake256_absorb(&state, mu, SEEDBYTES);
    shake256_absorb(&state, m, mlen);
    shake256_finalize(&state);
    shake256_squeeze(mu, CRHBYTES, &state);

    /* Matrix-vector multiplication; compute Az - c2^dt1 */
    poly_challenge(&cp, c);
    polyvecl_ntt(&z);
    polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
    poly_ntt(&cp);
    polyveck_shiftl(&t1);
    polyveck_ntt(&t1);
    polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);
    polyveck_sub(&w1, &w1, &t1);

    // define difference
    diff = (Q-mat[target_row].vec[target_col].coeffs[unpack_map[target_coeff]])%Q;

    // un-comment the following line to fine-tune performance (OMP_NESTED=TRUE might be necessary)
    //#pragma omp parallel for num_threads(4) schedule(static) // if this line is not commented out, please remove the break statement in this loop
    for (int candidate_ntty = 0; candidate_ntty < Q; candidate_ntty++)
    {
        if (foundflag) // slide out of parallelized loop
        {
            continue;
        }
        polyveck w1_tmp = w1; // copy to thread-local memory
        keccak_state local_state;
        uint8_t local_buf[K*POLYW1_PACKEDBYTES];
        uint8_t c2[SEEDBYTES];
        
        w1_tmp.vec[target_row].coeffs[unpack_map[target_coeff]] += ((((signed long)diff * (signed long)candidate_ntty) % Q) * 8265825L) % Q;
        w1_tmp.vec[target_row].coeffs[unpack_map[target_coeff]] %= Q;
        polyveck_reduce(&w1_tmp);
        polyveck_invntt_tomont(&w1_tmp);

        /* Reconstruct w1 */
        polyveck_caddq(&w1_tmp);
        polyveck_use_hint(&w1_tmp, &w1_tmp, &h);
        polyveck_pack_w1(local_buf, &w1_tmp);

        /* Call random oracle and verify challenge */
        shake256_init(&local_state);
        shake256_absorb(&local_state, mu, CRHBYTES);
        shake256_absorb(&local_state, local_buf, K*POLYW1_PACKEDBYTES);
        shake256_finalize(&local_state);
        shake256_squeeze(c2, SEEDBYTES, &local_state);
        
        if (memcmp(c, c2, SEEDBYTES) == 0)
        {
            #pragma omp critical
            foundflag = 1;
            recovered_ntty = candidate_ntty;
            break;
        }
    }
    res = recovered_ntty;
    if (recovered_ntty == 0) // find ineffective fault
    {
        #pragma omp critical
        cerr << "ineffective! ----------------------------------------------\n";
        return 0;
    }
    return foundflag;
}

int recover_s1_fA_ntt_coeff(unsigned int target_row,
                    unsigned int target_col, 
                    unsigned int target_coeff, 
                    const uint8_t *pk,
                    const polyvecl mat[K])
{
    uint8_t sm[10 + CRYPTO_BYTES];
    size_t smlen;
    ZZ_p ntty_coeff(0), nttz_coeff(0), nttc_coeff(0), ntts1_coeff(0);
    int ret;
    uint8_t c[SEEDBYTES];
    poly cp;
    polyvecl z;
    polyveck h;
    FILE *f;
    char fname[256];
    map<int, int> results;

    for (size_t cnt = 0; cnt < NFAULTS; cnt++) // iterate over faults at this position
    {
        // read signature
        snprintf(fname, 256, "%s/r%d_c%d/coeff%d/sig%lu", dirname, target_row, target_col, target_coeff, cnt);
        f = fopen(fname, "r");
        if (f == NULL)
        {
            #pragma omp critical
            cerr << "could not open: " << fname << "\n";
            continue;
        }
        fseek(f, 0, SEEK_END);
        fseek(f, -(int)CRYPTO_BYTES-10, SEEK_CUR);
        smlen = fread(sm, sizeof(uint8_t), 10 + CRYPTO_BYTES, f); // read signature and 10 bytes of message from end of file
        fclose(f);
        if (smlen <= CRYPTO_BYTES)
        {
            #pragma omp critical
            cerr << "too few data: " << fname << "\n";
            continue;
        }

        //get one coeff of the ntt rep. of y by correcting the signature
        if (!correct_sig_fA(ntty_coeff, sm + smlen - CRYPTO_BYTES - 10, CRYPTO_BYTES, sm + smlen - 10, 10, pk, mat, target_row, target_col, target_coeff, 1))
        {
            continue;
        }

        //now calculate the corresponding coeff in the ntt rep. of s1
        ret = unpack_sig(c, &z, &h, sm + smlen - CRYPTO_BYTES - 10);
        if (ret != 0)
        {
            #pragma omp critical
            cerr << "sig unpacking failed in recover_s1_fA_ntt_coeff\n";
            return -1;
        }
        polyvecl_ntt(&z);

        poly_challenge(&cp, c);
        poly_ntt(&cp);

        nttz_coeff = z.vec[target_col].coeffs[unpack_map[target_coeff]];
        nttc_coeff = cp.coeffs[unpack_map[target_coeff]];
        ntts1_coeff = nttz_coeff - ntty_coeff;

        div(ntts1_coeff, ntts1_coeff, nttc_coeff);
        conv(ret, ntts1_coeff);

        if (results.find(ret) == results.end()) // if this candidate is not yet contained in the map
        {
            results[ret] = 1;
        } else { // if it is already contained
            results[ret] += 1;
        }
        if (results[ret] >= NFAULTS_MINIMUM) // we have found our winner!
        {
            return ret;
        }
    }
    return -1; // return not found if NFAULTS_MINIMUM was never reached
}

void recover_s1_fA(const uint8_t *pk)
{
    polyvecl s1_recovered = {-1}, mat[K];
    polyveck t;
    uint8_t rho[SEEDBYTES];
    unpack_pk(rho, &t, pk);
    polyvec_matrix_expand(mat, rho);
    int total = 0, found = 0;

    //get all ntt(s1) coeff with recover_s1_fA_coeff
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (auto i = 0; i < L; i++)
    {
        for (auto j = 0; j < N; j++)
        {

            init_ntl(); // init within each thread!
            int tmp = recover_s1_fA_ntt_coeff(0, i, j, pk, mat);
            #pragma omp critical
            {
                if (tmp == -1)
                {
                    cout << "could not recover " << i << ", " << j;
                } else {
                    cout << "RECOVERED " << i << ", " << j << ": " << tmp;
                    s1_recovered.vec[i].coeffs[j] = tmp;
                    found += 1;
                }
                total += 1;
                cout << " (found " << found << " of " << total << ")\r" << flush;
            }
        }  
    }
    cout << "\n";
    char fname[256];
    snprintf(fname, 256, "%s/s1_recovered", dirname);
    FILE *f = fopen(fname, "w");
    if (f == NULL || fwrite(&s1_recovered, sizeof(polyvecl), 1, f) != 1)
    {
        cerr << "could not open or write " << fname << ", but here is the key I recovered: \n";
        for (auto i = 0; i < L; i++)
        {
            for (auto j = 0; j < N; j++)
            {
                cerr << s1_recovered.vec[i].coeffs[j] << ", ";
            }
            cerr << "\n";
        }
        return;
    }
    fclose(f);
}

int main(void)
{
    chrono::seconds elapsed, min = 0s, max = 0s;
    FILE *f;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    char fname[256];

    // iterate over signature sets (each set has its own key pair)
    for (auto mset = 0; mset < NSETS; mset++)
    {
        auto starttime = chrono::system_clock::now();

        snprintf(dirname, 128, "dilithium%d-sim%d/", DILITHIUM_MODE, mset); // change directory name HERE if necessary!

        // read pk
        snprintf(fname, 256, "%spk", dirname);
        f = fopen(fname, "r");
        if (f == NULL)
        {
            printf("file open error %s\n", fname);
            return -1;
        }
        if (fread(pk, 1, CRYPTO_PUBLICKEYBYTES, f) != CRYPTO_PUBLICKEYBYTES)
        {
            printf("bad pk len.\n");
            return -1;
        }
        fclose(f);
        // end read pk

        recover_s1_fA(pk); // call recovery code
        auto endtime = chrono::system_clock::now();
        auto elapsedtime = chrono::duration_cast<chrono::seconds>(endtime-starttime);

        // update min/max of elapsedtime
        if (mset == 0)
        {
            elapsed = elapsedtime;
            min = elapsedtime;
            max = elapsedtime;
        } else {
            elapsed += elapsedtime;
            if (elapsedtime > max)
            {
                max = elapsedtime;
            }
            if (elapsedtime < min)
            {
                min = elapsedtime;
            }
        }
        cout << elapsedtime.count() << "\n";
    }
    elapsed /= NSETS;
    cout << "avg: " << elapsed.count() << "\nmin: " << min.count() << "\nmax: " << max.count() << "\n";
    return 0;
}