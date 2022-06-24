# compilation flags
CC=/usr/bin/cc
CFLAGS=-O3 -Wall -std=c99 -g

CXX=/usr/bin/c++
CXX_FLAGS=-std=c++11 -Wall -Wextra -DNDEBUG
CXX_OPT_FLAGS=-O3 -ffast-math -funroll-loops -msse4.2 -march=native -DHAVE_CXA_DEMANGLE

LIBS=-lsdsl -ldivsufsort -ldivsufsort64

EXECS=pfwg.x newscanNT.x tfm_index_construct.x tfm_index_invert.x

.PHONY: all test clean
all: $(EXECS)

test: $(EXECS)
	./bigbwt -w 4 -p 50 -i data/yeast.raw

clean:
	rm -f *.o gsa/*.o
	rm -f *.x

# libs
gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

# executables
newscanNT.x: newscan.cpp utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

pfwg.x: pfwg.cpp tfm_index.hpp  gsa/gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ pfwg.cpp -Iinclude -Llib gsa/gsacak.o utils.o -pthread -ldl $(LIBS)

tfm_index_construct.x: tfm_index_construct.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_OPT_FLAGS) $(C_OPTIONS) \
	-Iinclude -Llib tfm_index_construct.cpp -o tfm_index_construct.x $(LIBS)

tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_OPT_FLAGS) $(C_OPTIONS) \
	-Iinclude -Llib tfm_index_invert.cpp -o tfm_index_invert.x $(LIBS)

num_runs: num_runs.c
	$(CC) num_runs.c -o num_runs

