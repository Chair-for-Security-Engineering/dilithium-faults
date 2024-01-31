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



test_randomstuff: test_latticered.cpp lattice_red.cpp lattice_red.h .dilithium2
	$(CXX) $(CXXFLAGS) $(CORR_BOUND) -DDILITHIUM_MODE=2 *.o lattice_red.cpp test_latticered.cpp -o test_latticered -lmpfr -lntl -lgmp `pkg-config --libs fplll`


.dilithium2: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) -std=c11 -c


clean:
	rm -f test_latticered
	rm -f *.o
