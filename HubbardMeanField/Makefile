CXX = mpicxx

INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp

LDLIBS = $(LINALG_FLAG)

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

HBBRD_SRCS=Model.cpp BasicHubbardModel.cpp HubbardCDW.cpp UsingBroyden.cpp
UTIL_SRCS=InputFileReader.cpp OutputWriter.cpp

PART_SRCS=
SRCS=$(addprefix Utility/, $(UTIL_SRCS)) $(addprefix Hubbard/, $(HBBRD_SRCS)) $(PART_SRCS) Hubbard_Mean_Field.cpp

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

	
