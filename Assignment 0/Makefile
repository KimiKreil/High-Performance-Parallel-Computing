CXX := g++
CXXFLAGS := -g -std=c++11

sir: sir.o
	$(CXX) sir.o -o sir # Runs second ... with object files you only recompile changes

sir.o: sir.cpp
	$(CXX) -c sir.cpp $(CXXFLAGS) -o sir.o # Runs first

clean:
	rm -fr sir.o sir sir_output.txt 
			
