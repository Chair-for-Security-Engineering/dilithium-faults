DILITHIUM_PATH := dilithium/avx2

LIBS := -lntl -lgmp -l:libpqcrystals_fips202_avx2.so -l:libpqcrystals_fips202x4_avx2.so
LIB_PATH := -L$(DILITHIUM_PATH)
CXXFLAGS := -O3 -Wall -fopenmp -march=native

all: test_2_fA test_3_fA test_5_fA gensig2 gensig3 gensig5

$(DILITHIUM_PATH):
	$(MAKE) -C $@ -j shared

test_2_fA: correction_attacks.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=2 -DFAULT=2 correction_attacks.cpp -o test_2_fA -l:libpqcrystals_dilithium2_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)
test_3_fA: correction_attacks.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=3 -DFAULT=2 correction_attacks.cpp -o test_3_fA -l:libpqcrystals_dilithium3_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)
test_5_fA: correction_attacks.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=5 -DFAULT=2 correction_attacks.cpp -o test_5_fA -l:libpqcrystals_dilithium5_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)

gensig2: gen_faulty_sigs.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=2 -DFAULT=2 gen_faulty_sigs.cpp -o gen_faulty_sig2 -l:libpqcrystals_dilithium2_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)
gensig3: gen_faulty_sigs.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=3 -DFAULT=2 gen_faulty_sigs.cpp -o gen_faulty_sig3 -l:libpqcrystals_dilithium3_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)
gensig5: gen_faulty_sigs.cpp $(DILITHIUM_PATH)
	$(CXX) $(CXXFLAGS) -DDILITHIUM_MODE=5 -DFAULT=2 gen_faulty_sigs.cpp -o gen_faulty_sig5 -l:libpqcrystals_dilithium5_avx2.so $(LIBS) $(INC_PATH) $(LIB_PATH)

clean:
	rm -rf *.o test*fA

.PHONY: all $(DILITHIUM_PATH)
