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


#include <cstring>
#include <fplll.h>


using namespace fplll;


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



/****************************/
/*ntl <-> fplll conversions */
/****************************/
void ntl_to_fplll_mat(ZZ_mat<mpz_t> &out, const Mat<ZZ_p> &in)
{
  out.set_cols(in.NumCols());
  out.set_rows(in.NumRows());
  for (auto i = 0; i < in.NumCols(); i++)
  {
    for (auto j = 0; j < in.NumRows(); j++)
    {
      out[j][i] = conv<long>(rep(in[j][i]));
    }
  }
}

void fplll_to_ntl_matrow(Vec<ZZ_p> &out, const MatrixRow<Z_NR<mpz_t>> &in)
{
  out.SetLength(in.size());
  for (auto i = 0; i < in.size(); i++)
  {
    out[i] = ZZ_p(in[i].get_si());
  }
}




/**************************************/
/* calculate the INTT Matrix          */
/* the rows form the lattice basis    */
/**************************************/
Mat<ZZ_p> get_intt_mat()
{
  Mat<ZZ_p> INTT_mat;
  INTT_mat.SetDims(N,N);
  poly s1_hat;

  for (auto i = 0; i < N; i++)
  {
      for (auto j = 0; j < N; j++)
      {
          s1_hat.coeffs[j] = 0;
      }
      s1_hat.coeffs[i] = 1;
        
      poly_invntt(&s1_hat);
      to_ntl(INTT_mat[i], s1_hat);
  }
  transpose(INTT_mat);
  return INTT_mat;
}



/*************************************************************/
/* make the reduce intt matrix for a single poly entry       */
/* in: Vec<ZZ_p> &known_coeff = known s_hat coeffs in order  */
/*     Vec<ZZ> &known_coeff_pos = 0 for unknown, 1 for known */
/*     INTT_mat                                              */
/*************************************************************/
void build_red_mat_single(Mat<ZZ_p> &INTT_red, const Vec<ZZ_p> &known_coeff, const Vec<int> &known_coeff_pos, const int nr_known, const Mat<ZZ_p> &INTT_mat)
{
  INTT_red.SetDims(N-nr_known+1,N);
  clear(INTT_red);
  int j = 0;
  int k = 0;
  for (auto i = 0; i < N; i++)
  {
    /*if the coeff is known, mult with corresponding INTT row and add to last row*/
    if (known_coeff_pos[i] == 1)
    {
      INTT_red[N-nr_known] += known_coeff[k]*INTT_mat[i];
      k++;
    }
    else if (known_coeff_pos[i] == 0)
    {
      /*if the coeff is unknown, just copy the INTT row*/
      INTT_red[j] = INTT_mat[i];
      j++;
    }
    else
    {
      throw logic_error("known_coeff_pos is not just 0, 1");
    }
  }
}




/************************************************/
/*make reduced echelon form of the given matrix */
/*without the last row                          */
/************************************************/
void make_red_echelon(Mat<ZZ_p> &A)
{
  long rows = A.NumRows();
  long colmns = A.NumCols();
  Mat<ZZ_p> A_wo_last;
  A_wo_last = A;
  A_wo_last.SetDims(rows-1, colmns);
  gauss(A_wo_last);
  if (rows > colmns)
  {
    cout << "make_red_echelon: are you sure you have transposed?" << endl;
  }
  //invert vorderen teil, mult drauf
  Mat<ZZ_p> inv_left_part, temp;
  inv_left_part.SetDims(rows-1, rows-1);
  temp = transpose(A_wo_last);
  temp.SetDims(rows-1, rows-1);
  temp = transpose(temp);
  inv(inv_left_part, temp);
  A_wo_last = inv_left_part * A_wo_last;
  for (auto i = 0; i < rows-1; i++)
  {
    A[i] = A_wo_last[i];
  }
  
}



/**************************************************/
/*embedding                                       */
/*input is supposed to be in red echelon row form */
/*aka left side is identity, last row is sideinfo */
/**************************************************/
void embedd(ZZ_mat<mpz_t> &A, const int nr_known)
{
  A.set_cols(N+1);
  for (auto i = 0; i < N-nr_known; i++)
  {
    A[i][N] = 0;
  }
  A[N-nr_known][N] = 1;

  A.set_rows(N+1);
  for (auto i = N-nr_known+1; i < N+1; i++)
  {
    for (auto j = 0; j < N+1; j++)
    {
      A[i][j] = 0;
    }
    //cout << i << endl;
    A[i][i-1] = Q;
  }
  
}


/**********************************/
/*main lattice_reduction function */
/*combines the functions above    */
/**********************************/
int lattice_reduction(Vec<ZZ_p> &found_sol, const Vec<ZZ_p> &known_coeff, const Vec<int> &known_coeff_pos, const int nr_known, int block_size)
{
  Mat<ZZ_p> INTT_red;
  Mat<ZZ_p> INTT_mat;
  INTT_mat = get_intt_mat();
  
  build_red_mat_single(INTT_red, known_coeff, known_coeff_pos, nr_known, INTT_mat);
  make_red_echelon(INTT_red);
  ZZ_mat<mpz_t> INTT_red_fplll;
  ntl_to_fplll_mat(INTT_red_fplll, INTT_red);
  embedd(INTT_red_fplll, nr_known);
  int status = 0;

  status = bkz_reduction(INTT_red_fplll, block_size, BKZ_NO_LLL, FT_MPFR, 128);
  fplll_to_ntl_matrow(found_sol, (INTT_red_fplll[0]));

  return status;
}



/*****************************************/
/*a function for testing the lattice red */
/*with random known shat coefficients    */
/*****************************************/
int test_latticered(int block_size, int nr_known)
{
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  crypto_sign_keypair(pk, sk);
  
  uint8_t rho[SEEDBYTES];
  uint8_t tr[SEEDBYTES];
  uint8_t key[SEEDBYTES];
  polyveck t0;
  polyvecl s1; 
  polyveck s2;
  Vec<ZZ_pX> s1_orig;
  //get the original key
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
  to_ntl(s1_orig, s1);

  //get some shat coeffs
  polyvecl_ntt(&s1);
  polyvecl_reduce(&s1);

  Vec<ZZ_p> known_coeff;
  known_coeff.SetLength(nr_known);
  
  Vec<int> known_coeff_pos;
  known_coeff_pos.SetLength(N);

  int success = 0;
  /*for all s1 components */
  for (auto comp = 0; comp < L; comp++)
  {  
    clear(known_coeff);
    for (auto i = 0; i < N; i++)
    {
      known_coeff_pos[i] = 0;
    }

    for (auto i = 0; i < nr_known; i++)
    {
      known_coeff_pos[i] = 1;
    }
    random_shuffle(known_coeff_pos.begin(), known_coeff_pos.end());
    int temp = 0;
    for (auto i = 0; i < N; i++)
    {
      if (known_coeff_pos[i] == 1)
      {
        known_coeff[temp] = s1.vec[comp].coeffs[i];
        temp++;
      }
      
    }


    Vec<ZZ_p> found_sol;
    found_sol.SetLength(N);
    clear(found_sol);

    auto tic = chrono::high_resolution_clock::now();
    try
    {
      success |= lattice_reduction(found_sol, known_coeff, known_coeff_pos, nr_known, block_size);
    }
    catch(const std::exception& e)
    {
      std::cerr << e.what() << '\n';
      success |= 1;
    }
    auto toc = chrono::high_resolution_clock::now();
    auto secs = sf(toc-tic).count();

    Vec<ZZ_p> correct_sol;
    correct_sol = s1_orig[comp].rep;
    correct_sol.append(ZZ_p(1));

    if (found_sol == correct_sol)
    {
      cout << "Found the correct s1 component in " << secs << " secs." << endl;
    }
    else
    {
      cout << "The found solution =/= correct solution after " << secs << " secs." << endl;
    }
  }
  return success;
}

