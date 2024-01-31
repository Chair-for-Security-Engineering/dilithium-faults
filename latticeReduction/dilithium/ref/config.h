#ifndef CONFIG_H
#define CONFIG_H

//#define DILITHIUM_MODE 2
//#define DILITHIUM_USE_AES
#define DILITHIUM_RANDOMIZED_SIGNING //all attacks are targeting the randomized variant
//#define USE_RDPMC
//#define DBENCH

/* NEW =======================================================================*/
#ifndef SET_Y_ZERO
#define SET_Y_ZERO 1 //corresponds to the original ref implementation
#endif

#ifndef SHUFFLING
#define SHUFFLING 0
#endif
/*============================================================================*/

#ifndef DILITHIUM_MODE
#define DILITHIUM_MODE 2
#endif

#ifdef DILITHIUM_USE_AES
#if DILITHIUM_MODE == 2
#define CRYPTO_ALGNAME "Dilithium2-AES"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium2aes_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium2aes_ref_##s
#elif DILITHIUM_MODE == 3
#define CRYPTO_ALGNAME "Dilithium3-AES"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium3aes_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium3aes_ref_##s
#elif DILITHIUM_MODE == 5
#define CRYPTO_ALGNAME "Dilithium5-AES"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium5aes_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium5aes_ref_##s
#endif
#else
#if DILITHIUM_MODE == 2
#define CRYPTO_ALGNAME "Dilithium2"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium2_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium2_ref_##s
#elif DILITHIUM_MODE == 3
#define CRYPTO_ALGNAME "Dilithium3"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium3_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium3_ref_##s
#elif DILITHIUM_MODE == 5
#define CRYPTO_ALGNAME "Dilithium5"
#define DILITHIUM_NAMESPACETOP pqcrystals_dilithium5_ref
#define DILITHIUM_NAMESPACE(s) pqcrystals_dilithium5_ref_##s
#endif
#endif

#endif
