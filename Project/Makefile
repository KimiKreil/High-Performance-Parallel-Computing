CXX := g++

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -ffast-math

# Add profiling to code
#OPT=-O1 -pg

# Faster compilation time
#OPT=-O1

CXXFLAGS := $(OPT) -Wall -march=znver1 -std=c++14

default: seq vec #sol flem

seq: Ising_Model_seq.cpp
	$(CXX) Ising_Model_seq.cpp $(CXXFLAGS) -o seq

vec: Ising_Model_vec.cpp
	$(CXX) Ising_Model_vec.cpp $(CXXFLAGS) -fopenmp -o vec


clean:
	rm -fr seq vec # vec sol flem
