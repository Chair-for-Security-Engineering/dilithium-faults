#!/bin/sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./dilithium/avx2
export OMP_NESTED=TRUE

./gen_faulty_sig2
./test_2_fA

./gen_faulty_sig3
./test_3_fA

./gen_faulty_sig5
./test_5_fA