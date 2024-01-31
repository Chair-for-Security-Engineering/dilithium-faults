#ifndef SIGN_H
#define SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

#define challenge DILITHIUM_NAMESPACE(challenge)
void challenge(poly *c, const uint8_t seed[SEEDBYTES]);

#define crypto_sign_keypair DILITHIUM_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature DILITHIUM_NAMESPACE(signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);
#define crypto_sign_signature_skip DILITHIUM_NAMESPACE(signature_skip)
int crypto_sign_signature_skip(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk,
                          const uint32_t target_comp,
                          const uint32_t target_coeff);
#define crypto_sign_signature_fA DILITHIUM_NAMESPACE(signature_fA)
int crypto_sign_signature_fA(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk,
                          const uint32_t target_row,
                          const uint32_t target_col,
                          const uint32_t target_coeff,
                          const uint32_t delta);

#define crypto_sign DILITHIUM_NAMESPACETOP
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);
#define crypto_sign_skip DILITHIUM_NAMESPACE(skip)
int crypto_sign_skip(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk,
                const uint32_t target_comp,
                const uint32_t target_coeff);
#define crypto_sign_fA DILITHIUM_NAMESPACE(fA)
int crypto_sign_fA(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk,
                const uint32_t target_row,
                const uint32_t target_col,
                const uint32_t target_coeff,
                const uint32_t delta);


#define crypto_sign_verify DILITHIUM_NAMESPACE(verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open DILITHIUM_NAMESPACE(open)
int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#endif
