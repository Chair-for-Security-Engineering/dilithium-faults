#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <unistd.h>
#include <sys/syscall.h>

#include "NTL/lzz_p.h"
#include "NTL/vec_lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
#include "NTL/mat_lzz_p.h"
#include "NTL/pair.h"

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
}

using namespace NTL;

//fault type: 1 for skipping
#ifndef FAULT
#define FAULT 1 
#endif

//number of tries before the brute-force correction is aborted
#ifndef CORR_TRIES
#define CORR_TRIES 2*BETA
#endif

#ifndef MLEN
#define MLEN 59
#endif

//global variables for statistical values
extern int faults;


/******************************************************************************************************/
/* helper functions                                                                                   */
/******************************************************************************************************/

void init_ntl();
void to_ntl(ZZ_pX& out, const poly& in);
void to_ntl(Vec<ZZ_p>& out, const poly& in);
void to_ntl(Vec<ZZ_pX>& out, const polyvecl& in);
void from_ntl(poly& out, const ZZ_pX& in);
void from_ntl(poly& out, const Vec<ZZ_p>& in);
void from_ntl(polyvecl& out, const Vec<ZZ_pX>& in);

void eq_from_sig(Vec<ZZ_p> *eq_coeffs, ZZ_p *eq_sol, unsigned int *target_comp_actual, uint8_t *sm, uint8_t *m, unsigned int target_comp, unsigned int target_coeff, const uint8_t *pk, size_t mlen);
void coeffs_weird_rot(Vec<ZZ_p>& in, long rot);
Vec<Vec<ZZ_p>> get_indices(int n, int t, int d, bool max_incl);



/******************************************************************************************************/
/* skipping fault                                                                                     */
/******************************************************************************************************/

Vec<ZZ_pX> recover_s1_skip(const uint8_t *sk, const uint8_t *pk, size_t mlen);
void recover_s1_skip_comp(Vec<ZZ_p> *solution, const unsigned int target_component, const uint8_t *sk, const uint8_t *pk, size_t mlen);
int correct_sig_skip(uint8_t *sm, 
    	            uint8_t *m,
                    unsigned int target_comp, 
                    unsigned int target_coeff, 
                    const uint8_t *pk,
                    size_t mlen);



/******************************************************************************************************/
/* shuffled skipping fault                                                                            */
/******************************************************************************************************/

Vec<ZZ_pX> recover_s1_skip_shuff(const uint8_t *sk, const uint8_t *pk, size_t mlen);
int correct_sig_skip_shuff(unsigned int *target_comp_actual,
        	    	unsigned int *target_coeff_actual,
                    uint8_t *sm, 
                    uint8_t *m,
                    const uint8_t *pk,
                    size_t mlen);



/******************************************************************************************************/
/* forgery                                                                                            */
/******************************************************************************************************/

int crypto_forge_signature(uint8_t *sig, size_t *siglen, const uint8_t *m, size_t mlen, const Vec<ZZ_pX> s1_rec, const uint8_t *pk);
int crypto_forge_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen, const Vec<ZZ_pX> s1_rec, const uint8_t *pk);
