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

#include <fplll.h>

using namespace NTL;
using namespace fplll;


void init_ntl();
void to_ntl(ZZ_pX& out, const poly& in);
void to_ntl(Vec<ZZ_p>& out, const poly& in);
void from_ntl(poly& out, const ZZ_pX& in);
void from_ntl(poly& out, const Vec<ZZ_p>& in);
void to_ntl(Vec<ZZ_pX>& out, const polyvecl& in);
void from_ntl(polyvecl& out, const Vec<ZZ_pX>& in);
Mat<ZZ_p> get_intt_mat();
void ntl_to_fplll_mat(ZZ_mat<mpz_t> &out, Mat<ZZ_p> &in);
void make_red_echelon(Mat<ZZ_p> &A);
void build_red_mat_single(Mat<ZZ_p> &INTT_red, 
                            const Vec<ZZ_p> &known_coeff, 
                            const Vec<ZZ> &known_coeff_pos, 
                            const int nr_known,
                            const Mat<ZZ_p> &INTT_mat);
void embedd(ZZ_mat<mpz_t> &A, const int nr_known);
int lattice_reduction(Vec<ZZ_p> &found_sol, const Vec<ZZ_p> &known_coeff, const Vec<int> &known_coeff_pos, const int nr_known, int block_size);
int test_latticered(int block_size, int nr_known);