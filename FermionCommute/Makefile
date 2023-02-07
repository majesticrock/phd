CXX = g++

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp -DNDEBUG

LDLIBS = $(LINALG_FLAG)

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

PART_SRCS=
SRCS=$(PART_SRCS) FermionCommute.cpp

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main 

build/main: $(OBJS) | build
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build
	mkdir -p build/Hubbard
	mkdir -p build/Utility

clean:
	rm -f $(OBJS)
	rm -f main
	rm -rf build

	
