#include "correction_attacks.h"

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
#include "dilithium/ref/poly.h"
#include "dilithium/ref/polyvec.h"
#include "dilithium/ref/packing.h"
#include "dilithium/ref/sign.h"
#include "dilithium/ref/randombytes.h"
#include "dilithium/ref/params.h"
#include "dilithium/ref/fips202.h"
}


using namespace std;
using namespace NTL;

typedef chrono::duration<double, ratio< 1>>   sf;



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
/* coeffs_weird_rot anti-cyclic rotates a given vector                         */
/* used to get coeff of a polynom*polynom multiplication result                */
/* returns the coeffs in the corresponding row of the matrix (not column)      */
/*                                                                             */
/* Only used for the skipping fault attack                                     */
/*******************************************************************************/

void coeffs_weird_rot(Vec<ZZ_p>& in, long rot)
{
    Vec<ZZ_p> temp;
    temp.SetLength(in.length());
    //reverse the input vector
    for (auto i = 0; i < in.length(); i++)
    {
        temp[i] = in[N-1-i];
    }
    
    //rotate right rot+1 times and -1 the entries that did not "cycle back"
    rot += 1;
    for (auto i = 0; i < in.length()-rot; i++)
    {
        in[i+rot] = (-1)*temp[i];
    }
    for (auto i = 0; i < rot; i++)
    {
        in[i] = temp[N-rot+i];
    }
}



/*******************************************************************************/
/* eq_from_sig extracts an equation from a signature                           */
/* takes a signed message, a target component and a target coefficient         */
/* and returns the underlying equation in the corresponding s1 component       */
/*                                                                             */
/* respects SET_Y_ZERO                                                         */
/* respects SHUFFLING in the SET_Y_ZERO = 0 case                               */
/*                                                                             */
/* Only used for the skipping fault attack                                     */
/*******************************************************************************/

void eq_from_sig(Vec<ZZ_p> *eq_coeffs, 
                ZZ_p *eq_sol,
                unsigned int *target_comp_actual,
                uint8_t *sm, 
                uint8_t *m,
                unsigned int target_comp, 
                unsigned int target_coeff,
                const uint8_t *pk,
                size_t mlen)
{
    int ret;
    uint8_t c[SEEDBYTES];
    polyvecl z;
    polyveck h; 
    ZZ_pX temp_z;
    poly temp_coeffs;
    unsigned int target_coeff_actual;
    target_coeff_actual = target_coeff;
    *target_comp_actual = target_comp;
    
    //get the signature components, especially c and z
    ret = unpack_sig(c, &z, &h, sm);
    if (ret != 0)
    {
        throw logic_error("sig unpacking failed in eq_from_sig");
        return;
    }
    
    //the cases of SET_Y_ZERO only differ in the equation solution
    if (SET_Y_ZERO == 1) //z=cs1 case
    {
        //if the addition of y is skipped, z gives the equation solution directly
        to_ntl(temp_z, z.vec[*target_comp_actual]);
        *eq_sol = coeff(temp_z, target_coeff_actual);
    }
    else //z=y case
    {
        //find cs1 by correcting the signature
        //the correction value is the equation solution
        if (SHUFFLING == 0)
        {
            try
            {
                *eq_sol = correct_sig_skip(sm, m, *target_comp_actual, target_coeff_actual, pk, mlen);
            }
            catch(const exception& e)
            {
                throw logic_error(e.what());
            }
        }
        else //shuffling case
        {
            try
            {
                *eq_sol = correct_sig_skip_shuff(target_comp_actual, &target_coeff_actual, sm, m, pk, mlen);
            }
            catch(const exception& e)
            {
                throw logic_error(e.what());
            }
        }
    }

    //c gives the equation coefficients after rotating them with the anti-cyclic rotation
    poly_challenge(&temp_coeffs, c);
    to_ntl(*eq_coeffs, temp_coeffs);
    coeffs_weird_rot(*eq_coeffs, target_coeff_actual); 
}




/*******************************************************************************/
/* SKIPPING FAULT                                                              */
/*                                                                             */
/* correct_sig_skip takes a skip-faulted signature and the public key          */
/* then the faulted coeff is corrected to yield a valid signature              */
/* the difference (correction value) is given as the output (int)              */
/*                                                                             */
/* recover_s1_skip_comp recovers one component (polynom) of s1                 */
/* it needs access to a signing oracle crypto_sign_skip(sk)                    */
/* and returns the recovered polynom as vector of its coefficients             */
/*                                                                             */
/* recover_s1_skip simply calls recover_s1_skip_comp for all L components      */
/* to get the full s1                                                          */
/*******************************************************************************/

int correct_sig_skip(uint8_t *sm, 
                    uint8_t *m,
                    unsigned int target_comp, 
                    unsigned int target_coeff, 
                    const uint8_t *pk,
                    size_t mlen)
{
    int corrected_by;
    corrected_by = 0;
    
    //unpack the signature
    int ret;
    uint8_t c[SEEDBYTES];
    polyvecl z;
    polyveck h; 

    ret = unpack_sig(c, &z, &h, sm);
    if (ret != 0)
    {
        throw logic_error("sig unpacking failed in correct_sig_skip");
    }	

    int i, sign;
    i = 0;
    sign = -1; 

    //counter for how many values have been tested
    int count;
    count =0;

    //while verify == -1 add to z_comp_coeff until verify == 0
    //the program tries 0, 1, -1, 2, -2, 3, -3 and so on
    while (crypto_sign_verify(sm, CRYPTO_BYTES, m, mlen, pk) == -1)
    {
        i++;
        sign *= -1;
        corrected_by += sign*i;
        
        (z.vec[target_comp]).coeffs[target_coeff] += sign*i;
        pack_sig(sm, c, &z, &h);

        count++;
        if (count > CORR_TRIES)
        {
            throw logic_error("correction exceeded correction bound in correct_sig_skip");
        }
        
    }
    return corrected_by;
}

void recover_s1_skip_comp(Vec<ZZ_p> *solution, const unsigned int target_component, const uint8_t *sk, const uint8_t *pk, size_t mlen)
{
    Mat<ZZ_p> eq_sys;
    eq_sys.SetDims(N,N);
    Vec<ZZ_p> eq_sys_line;

    Vec<ZZ_p> eq_sys_b;
    eq_sys_b.SetLength(N);
    ZZ_p eq_sys_b_comp;

    uint8_t m[mlen];
    uint8_t sm[mlen + CRYPTO_BYTES];
    size_t smlen;
    ZZ_p d;
    unsigned int target_comp_actual;
    target_comp_actual = target_component;

    //Generate N equations by using N random messages.
    new_eq_sys:
    for (auto i = 0; i < N; i++)
    {
        try_again:
        randombytes(m, mlen);

        //sign the message with the skipping fault, 
        //add the resulting linear equation to the system
        crypto_sign_skip(sm, &smlen, m, mlen, sk, target_comp_actual, i%N);
        faults++;
        try
        {
            eq_from_sig(&eq_sys_line, &eq_sys_b_comp, &target_comp_actual, sm, m, target_component, i%N, pk, mlen);
        }
        catch(const std::exception& e)
        {
            goto try_again;
        }
        eq_sys[i] = eq_sys_line;
        eq_sys_b[i] = eq_sys_b_comp;
    }

    //now solve the resulting linear equation system in the coefficients of the targeted s1 component
    if (CORR_TRIES == 0)
    {
        /*this case only generates equations with solution 0, the system to solve is Cx = 0
        hence this case needs to find the kernel of the system matrix */
        Mat<ZZ_p> eq_sys_t, ker;
        eq_sys_t = transpose(eq_sys);
        try
        {
            kernel(ker, eq_sys_t);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            cout << "kernel failed" << endl;
        }
        
        if (ker.NumRows() != 1)
        {
            throw logic_error("kernel had more or less than one dimension in recover_s1_skip_comp");
        }
        *solution = ker[0];
    }
    else //the case with non-zero solutions
    {
        solve(d, eq_sys, *solution, eq_sys_b);
        if (d == 0) //just in case the system was not solvable
        {
            //this is highly inefficient, because all eqations are thrown away
            //but it never happened during testing
            cout << "an equation system was not solvable in recover_s1_skip_comp" << endl;
            goto new_eq_sys;
        }
    }
}

Vec<ZZ_pX> recover_s1_skip(const uint8_t *sk, const uint8_t *pk, size_t mlen) 
{
    Vec<ZZ_p> solution;
    Vec<ZZ_pX> s1_recovered;
    solution.SetLength(N);
    s1_recovered.SetLength(L);

    if (CORR_TRIES != 0)
    {
        for (auto j = 0; j < L; j++)
        {
            int safety_stop = 0;
            try_again_1:
            try
            {
                recover_s1_skip_comp(&solution, j, sk, pk, mlen);
            }
            catch(const std::exception& e)
            {
                cout << e.what() << endl;
                safety_stop++;
                if (safety_stop > 10)
                {
                    throw logic_error("Something went wrong. s1 component recovery failed 10 times (recover_s1_skip).");
                }
                goto try_again_1;
            }
            
            for (auto i = 0; i < N; i++)
            {
                SetCoeff(s1_recovered[j], i, solution[i]);
            }
        }
    }
    else if (CORR_TRIES == 0)  //ineffective attack
    {
        Vec<Vec<int>> possible_factors;
        possible_factors.SetLength(L);

        for (auto j = 0; j < L; j++)
        {
            int safety_stop = 0;
            try_again_2:
            try
            {
                recover_s1_skip_comp(&solution, j, sk, pk, mlen);
            }
            catch(const std::exception& e)
            {
                cout << e.what() << endl;
                safety_stop++;
                if (safety_stop > 10)
                {
                    throw logic_error("Something went wrong. s1 component recovery failed 10 times (recover_s1_skip).");
                }
                
                goto try_again_2;
            }

            //filter the kernel factors
            int curr_length;
            possible_factors[j].SetLength(0);
            for (auto factor = 1; factor < Q/2 +1; factor++)
            {
                curr_length = possible_factors[j].length();
                for (auto i = 0; i < N; i++)
                {
                    if ((rep(solution[i]*ZZ_p(factor)) > ETA && rep(solution[i]*ZZ_p(factor))-Q < -ETA))
                    {
                        goto next_factor;
                    }        
                }
                possible_factors[j].SetLength(curr_length+2);
                possible_factors[j][curr_length] = factor;
                possible_factors[j][curr_length+1] = -factor;
                next_factor:
                continue;
            }

            if (possible_factors[j].length() < 1)
            {
                throw logic_error("kernel solution could not be turned into a valid s1 solution.");
            }

            for (auto i = 0; i < N; i++)
            {
                SetCoeff(s1_recovered[j], i, solution[i]);
            }
        }
  
        //Create an iterator of all remaining factor combinations 
        vector<vector<int>> iter = {{}};
        for (const auto& u : possible_factors) {
            vector<vector<int>> r;
            for (const auto& x : iter) {
                for (const auto y : u) {
                    r.push_back(x);
                    r.back().push_back(y);
                }
            }
            iter = move(r);
        }
            
        //Test the remaining factors by trying to forge a signature
        Vec<ZZ_pX> s1_recovered_temp = s1_recovered;
        for (const auto& fact_comb : iter)
        {
            for (auto comp = 0; comp < L; comp++)
            {
                s1_recovered_temp[comp] = s1_recovered[comp]*fact_comb[comp];
            }
            
            uint8_t m[MLEN];
            uint8_t sm[MLEN + CRYPTO_BYTES];
            size_t smlen;
            randombytes(m, MLEN);

            try
            {
                crypto_forge_sign(sm, &smlen, m, MLEN, s1_recovered_temp, pk);
            }
            catch(const std::exception& e)
            {
                continue;
            }
            return s1_recovered_temp;
        }
    }
    return s1_recovered;
}




/*******************************************************************************/
/* SKIPPING FAULT WITH SHUFFLING                                               */
/*                                                                             */
/* correct_sig_skip_shuff takes a skip-faulted signature and the public key    */
/* but has no knowledge of the targeted coeff.                                 */
/* It returns the correction value as well as the actually faulted coeff+comp  */
/*                                                                             */
/* recover_s1_skip_shuff attacks the scenario                                  */
/* in which the target comp+coeff is unknown                                   */
/* In this case, the comp.-wise brute-force is not a good idea (large overh.)  */
/* the brute-force is rather done by trying one value for all components       */
/* before trying the next                                                      */
/*******************************************************************************/

int correct_sig_skip_shuff(unsigned int *target_comp_actual,
        	    	unsigned int *target_coeff_actual,
                    uint8_t *sm, 
                    uint8_t *m,
                    const uint8_t *pk,
                    size_t mlen)
{
    int corrected_by;
    corrected_by = 0;
    
    //unpack the signature
    int ret;
    uint8_t c[SEEDBYTES];
    polyvecl z;
    polyveck h; 

    ret = unpack_sig(c, &z, &h, sm);
    if (ret != 0)
    {
        throw logic_error("sig unpacking failed in correct_sig_skip_shuff");
    }	

    int i, sign;
    i = 0;
    sign = -1; 

    //counter for how many values have been tested
    int count;
    count =0;
    int32_t temp;

    //while verify == -1 add to z_comp_coeff until verify == 0
    //the program tries 0, 1, -1, 2, -2, 3, -3 and so on
    while (crypto_sign_verify(sm, CRYPTO_BYTES, m, mlen, pk) == -1)
    {
        i++;
        sign *= -1;
        corrected_by += sign*i;
        
        //the value to be bruteforced is small,
        //so the program tries 1 for all possible coefficients, then -1 for all coefficients and so on
        for (auto target_comp_assumed = 0; target_comp_assumed < L; target_comp_assumed++)
        {
            for (auto target_coeff_assumed = 0; target_coeff_assumed < N; target_coeff_assumed++)
            {
                temp = (z.vec[target_comp_assumed]).coeffs[target_coeff_assumed];
                (z.vec[target_comp_assumed]).coeffs[target_coeff_assumed] += corrected_by;
                pack_sig(sm, c, &z, &h);
                if (crypto_sign_verify(sm, CRYPTO_BYTES, m, mlen, pk) == 0)
                {
                    *target_comp_actual = target_comp_assumed;
                    *target_coeff_actual = target_coeff_assumed;
                    return corrected_by;
                }
                (z.vec[target_comp_assumed]).coeffs[target_coeff_assumed] = temp;
            }        
        }

        count++;
        if (count > CORR_TRIES)
        {
            throw logic_error("correction exceeded correction bound in correct_sig_skip_shuff");
        }
        
    }
    return corrected_by;
}

Vec<ZZ_pX> recover_s1_skip_shuff(const uint8_t *sk, const uint8_t *pk, size_t mlen)
{
    Vec < ZZ_pX > s1_recovered;
    s1_recovered.SetLength(L);

    Vec < Mat<ZZ_p> > eq_sys_s;
    eq_sys_s.SetLength(L);
    Vec < Vec<ZZ_p> > eq_sys_b_s;
    eq_sys_b_s.SetLength(L);
    Vec < unsigned int > eq_counts;
    eq_counts.SetLength(L);
    for (auto i = 0; i < L; i++)
    {
        eq_sys_s[i].SetDims(N,N);
        eq_sys_b_s[i].SetLength(N);
        eq_counts[i] = 0;
    }
    
    Vec<ZZ_p> eq_sys_line;
    ZZ_p eq_sys_b_comp;

    uint8_t m[mlen];
    uint8_t sm[mlen + CRYPTO_BYTES];
    size_t smlen;
    ZZ_p d;
    unsigned int target_comp_actual;
    
    //Generate equations by using random messages until there are at least N equations for every component.
    new_eq_sys:
    for (auto i = 0; i < 2*L*N; i++)
    {
        //print progress
        if (i % 200 == 0)
        {
            cout << "generate equation " << i << endl;
        }
        
        try_again:
        randombytes(m, mlen);

        //sign the message with the skipping fault, 
        //add the resulting linear equation to the system for the recovered actual target component
        crypto_sign_skip(sm, &smlen, m, mlen, sk, 1, 1);
        faults++;
        try
        {
            eq_from_sig(&eq_sys_line, &eq_sys_b_comp, &target_comp_actual, sm, m, 1, 1, pk, mlen);
        }
        catch(const std::exception& e)
        {
            goto try_again;
        }

        if ((eq_sys_b_comp == 0) || (eq_counts[target_comp_actual] == N))
        {
            continue;
        }
        
        eq_sys_s[target_comp_actual][eq_counts[target_comp_actual]] = eq_sys_line;
        eq_sys_b_s[target_comp_actual][eq_counts[target_comp_actual]] = eq_sys_b_comp;
        eq_counts[target_comp_actual] += 1;

        int check_complete_systems = 0;
        for (auto comp = 0; comp < L; comp++)
        {
            if (eq_counts[comp] < N)
            {
                break;
            }
            else
            {
                check_complete_systems++;
            }
        }
        if (check_complete_systems == L)
            {
                break;
            }
    }

    //solve the resulting linear equation systems in the coefficients of the s1 components
    for (auto i = 0; i < L; i++)
    {
        Vec<ZZ_p> solution;
        solve(d, eq_sys_s[i], solution, eq_sys_b_s[i]);

        if (d == 0) //just in case the system was not solvable
        {
            //this is highly inefficient, because all eqations are thrown away
            //but it never happened in practice
            goto new_eq_sys;
        }

        for (auto j = 0; j < N; j++)
        {
            SetCoeff(s1_recovered[i], j, solution[j]);
        }   
    } 
    return s1_recovered;
}



/*******************************************************************************/
/* SIGNATURE FORGERY                                                           */
/*                                                                             */
/* Given the public key, s1 and a message, a verifiable signature is produced. */
/* Uses the original sign function definition as template.                     */
/*                                                                             */
/* For completeness the full sign proc. with the forgery is also implemented.  */
/*******************************************************************************/

int crypto_forge_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const Vec<ZZ_pX> s1_rec,
                          const uint8_t *pk)
{
    unsigned int n;
    uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
    uint8_t *rho, *tr, *key, *mu, *rhoprime;
    uint16_t nonce = 0;
    polyvecl mat[K], s1, y, z;
    polyveck t1, w1, w0, h, u;
    poly cp;
    keccak_state state;

    rho = seedbuf;
    tr = rho + SEEDBYTES;
    key = tr + SEEDBYTES;
    mu = key + SEEDBYTES;
    rhoprime = mu + CRHBYTES;
    unpack_pk(rho, &t1, pk);
    shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    
    /* Compute CRH(tr, msg) */
    shake256_init(&state);
    shake256_absorb(&state, tr, SEEDBYTES);
    shake256_absorb(&state, m, mlen);
    shake256_finalize(&state);
    shake256_squeeze(mu, CRHBYTES, &state);

    //Always use the randomized y generation. We dont know 'key' anyways
    randombytes(rhoprime, CRHBYTES);

    //Expand matrix and transform vectors
    polyvec_matrix_expand(mat, rho);
    from_ntl(s1, s1_rec);
    polyvecl_ntt(&s1);

    //Calculate u = As1 - t1*2^d = t0 - s2
    polyvec_matrix_pointwise_montgomery(&u, mat, &s1);
    polyveck_invntt_tomont(&u);
    polyveck_shiftl(&t1);
    polyveck_sub(&u, &u, &t1);
    polyveck_reduce(&u);

    polyveck_ntt(&u);

    int tries = 0;

    rej:
    tries++;
    
    /*safety stop*/
    if (tries > 500)
    {
        throw logic_error("forgery failed >500 times");
    }

    /* Sample intermediate vector y */
    polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

    /* Matrix-vector multiplication */
    z = y;
    polyvecl_ntt(&z);
    polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
    polyveck_reduce(&w1);
    polyveck_invntt_tomont(&w1);

    /* Decompose w and call the random oracle */
    polyveck_caddq(&w1);
    polyveck_decompose(&w1, &w0, &w1);
    polyveck_pack_w1(sig, &w1);

    shake256_init(&state);
    shake256_absorb(&state, mu, CRHBYTES);
    shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
    shake256_finalize(&state);
    shake256_squeeze(sig, SEEDBYTES, &state);
    poly_challenge(&cp, sig);
    poly_ntt(&cp);

    /* Compute z, reject if it reveals secret */
    polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
    polyvecl_invntt_tomont(&z);
    polyvecl_add(&z, &z, &y);
    polyvecl_reduce(&z);
    if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    {
        goto rej;
    }

    /* Compute hints for w1 using u */
    polyveck_pointwise_poly_montgomery(&h, &cp, &u);
    polyveck_invntt_tomont(&h);
    polyveck_reduce(&h);

    polyveck_add(&w0, &w0, &h);
    n = polyveck_make_hint(&h, &w0, &w1);
    if(n > OMEGA)
    {
        goto rej;
    }

    /* Write signature */
    pack_sig(sig, sig, &z, &h);
    *siglen = CRYPTO_BYTES;

    //Check, if the forged signature is valid
    if(crypto_sign_verify(sig, *siglen, m, mlen, pk) != 0)
    {
	    goto rej;
    }

    return tries;
}

int crypto_forge_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const Vec<ZZ_pX> s1_rec,
                const uint8_t *pk)
{
  size_t i;
  int tries;

  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  tries = crypto_forge_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, s1_rec, pk);
  *smlen += mlen;
  
  return tries;
}
