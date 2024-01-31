CC = gcc
CFLAGS = -Wall -Wextra -march=native -mtune=native -O3 -fomit-frame-pointer
NISTFLAGS = -O3 -fomit-frame-pointer -march=native -fPIC
SRCFOLDER = dilithium/ref
SOURCES = sign.c polyvec.c packing.c poly.c reduce.c ntt.c rounding.c fips202.c randombytes.c symmetric-shake.c
HEADERS = params.h sign.h polyvec.h packing.h poly.h reduce.h ntt.h rounding.h fips202.h randombytes.h symmetric.h
  
SOURCES := $(addprefix $(SRCFOLDER)/, $(SOURCES))
HEADERS := $(addprefix $(SRCFOLDER)/, $(HEADERS))

SOURCES2 = correction_attacks.cpp
HEADERS2 = correction_attacks.h

CXX = g++
CXXFLAGS = $(CFLAGS) -std=c++0x

###################################
# Change here for testing

CORR_BOUND = -DCORR_TRIES=2*BETA   #default 2*BETA
###################################


all: test_skip test_skip_shuff test_skip_ineff


.dilithium2: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) -std=c11 -c
.dilithium2_corr: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium2_corr_ineff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium2_corr_shuff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DSET_Y_ZERO=0 -DSHUFFLING=1 $(SOURCES) -std=c11 -c


.dilithium3: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 $(SOURCES) -std=c11 -c
.dilithium3_corr: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium3_corr_ineff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium3_corr_shuff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 -DSET_Y_ZERO=0 -DSHUFFLING=1 $(SOURCES) -std=c11 -c


.dilithium5: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 $(SOURCES) -std=c11 -c
.dilithium5_corr: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium5_corr_ineff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 -DSET_Y_ZERO=0 $(SOURCES) -std=c11 -c
.dilithium5_corr_shuff: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 -DSET_Y_ZERO=0 -DSHUFFLING=1 $(SOURCES) -std=c11 -c




#correction attack I using the skipping fault
test_skip: test_2_skip test_3_skip test_5_skip

test_2_skip: statistics_attack.cpp correction_attacks.cpp .dilithium2_corr
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=10 -DDILITHIUM_MODE=2 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_2_skip -lntl -lgmp
test_3_skip: statistics_attack.cpp correction_attacks.cpp .dilithium3_corr
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=10 -DDILITHIUM_MODE=3 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_3_skip -lntl -lgmp
test_5_skip: statistics_attack.cpp correction_attacks.cpp .dilithium5_corr
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=10 -DDILITHIUM_MODE=5 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_5_skip -lntl -lgmp



#correction attack I using the skipping fault, with shuffling
test_skip_shuff: test_2_skip_shuff test_3_skip_shuff test_5_skip_shuff

test_2_skip_shuff: statistics_attack.cpp correction_attacks.cpp .dilithium2_corr_shuff
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=2 -DDILITHIUM_MODE=2 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=1 *.o correction_attacks.cpp statistics_attack.cpp -o test_2_skip_shuff -lntl -lgmp
test_3_skip_shuff: statistics_attack.cpp correction_attacks.cpp .dilithium3_corr_shuff
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=2 -DDILITHIUM_MODE=3 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=1 *.o correction_attacks.cpp statistics_attack.cpp -o test_3_skip_shuff -lntl -lgmp
test_5_skip_shuff: statistics_attack.cpp correction_attacks.cpp .dilithium5_corr_shuff
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DTESTS=2 -DDILITHIUM_MODE=5 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=1 *.o correction_attacks.cpp statistics_attack.cpp -o test_5_skip_shuff -lntl -lgmp



#correction attack I using the skipping fault, with ineffective faults
test_skip_ineff: test_2_skip_ineff test_3_skip_ineff test_5_skip_ineff

test_2_skip_ineff: statistics_attack.cpp correction_attacks.cpp .dilithium2_corr_ineff
	$(CXX) $(CXXFLAGS) -DCORR_TRIES=0 -DTESTS=5 -DDILITHIUM_MODE=2 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_2_skip_ineff -lntl -lgmp
test_3_skip_ineff: statistics_attack.cpp correction_attacks.cpp .dilithium3_corr_ineff
	$(CXX) $(CXXFLAGS) -DCORR_TRIES=0 -DTESTS=5 -DDILITHIUM_MODE=3 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_3_skip_ineff -lntl -lgmp
test_5_skip_ineff: statistics_attack.cpp correction_attacks.cpp .dilithium5_corr_ineff
	$(CXX) $(CXXFLAGS) -DCORR_TRIES=0 -DTESTS=5 -DDILITHIUM_MODE=5 -DFAULT=1 -DSET_Y_ZERO=0 -DSHUFFLING=0 *.o correction_attacks.cpp statistics_attack.cpp -o test_5_skip_ineff -lntl -lgmp

	



clean:
	rm -f test_2_skip
	rm -f test_3_skip
	rm -f test_5_skip
	rm -f test_2_skip_shuff
	rm -f test_3_skip_shuff
	rm -f test_5_skip_shuff
	rm -f test_2_skip_ineff
	rm -f test_3_skip_ineff
	rm -f test_5_skip_ineff

	rm -f *.o
