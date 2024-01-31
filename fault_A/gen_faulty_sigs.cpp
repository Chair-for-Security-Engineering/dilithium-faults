#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cstdlib>
#include <time.h>

namespace fs = std::filesystem;

extern "C" {
#include "dilithium/avx2/poly.h"
#include "dilithium/avx2/polyvec.h"
#include "dilithium/avx2/packing.h"
#include "dilithium/avx2/sign.h"
#include "dilithium/avx2/randombytes.h"
#include "dilithium/avx2/params.h"
#include "dilithium/avx2/fips202.h"
}

using namespace std;

#define NFAULTS 10

void randombytes(uint8_t *out, size_t outlen)
{
  for (size_t i = 0; i < outlen; i++)
  {
    out[i] = rand();
  }
}

int main(void)
{
  uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES], msg[10] = {0}, sig[CRYPTO_BYTES];
  size_t siglen;
  ofstream f;
  char fname[256], dirname[256];
  srand(time(NULL));

  for (auto i = 0; i < NFAULTS; i++)
  {
    snprintf(dirname, 256, "dilithium%d-sim%d/", DILITHIUM_MODE, i);
    if ((!fs::is_directory(dirname) || !fs::exists(dirname)) && !fs::create_directory(dirname))
    {
      cerr << "Could not create directory " << dirname << "\n";
      continue;
    }

    crypto_sign_keypair(pk, sk);
    snprintf(fname, 256, "%spk", dirname);
    f.open(fname, ios::out | ios::binary);
    f.write((char*)pk, CRYPTO_PUBLICKEYBYTES);
    f.close();

    snprintf(fname, 256, "%ssk", dirname);
    f.open(fname, ios::out | ios::binary);
    f.write((char*)sk, CRYPTO_SECRETKEYBYTES);
    f.close();

    for (auto j = 0; j < L; j++)
    {
      snprintf(dirname, 256, "dilithium%d-sim%d/r0_c%d/", DILITHIUM_MODE, i, j);
      if ((!fs::is_directory(dirname) || !fs::exists(dirname)) && !fs::create_directory(dirname))
      {
        cerr << "Could not create directory " << dirname << "\n";
        continue;
      }
      for (auto k = 0; k < N; k++)
      {
        snprintf(dirname, 256, "dilithium%d-sim%d/r0_c%d/coeff%d/", DILITHIUM_MODE, i, j, k);
        if ((!fs::is_directory(dirname) || !fs::exists(dirname)) && !fs::create_directory(dirname))
        {
          cerr << "Could not create directory " << dirname << "\n";
          continue;
        }
        crypto_sign_signature_faulty(sig, &siglen, msg, 10, sk, j, k);
        snprintf(fname, 256, "%ssig0", dirname);
        f.open(fname, ios::out | ios::binary);
        f.write((char*)sig, CRYPTO_BYTES);
        f.write((char*)msg, 10);
        f.close();
      }
    }
  }
}