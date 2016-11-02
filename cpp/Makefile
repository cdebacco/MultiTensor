#CXXFLAGS=-Wall -O0 -g -Wno-deprecated -lgomp -fopenmp
CXXFLAGS=-Wall -std=c++0x -O3 -g -Wno-deprecated
LIBS=-L/usr/local/Cellar/boost/1.58.0/lib/ -lboost_program_options
LIBS2=-I /usr/include
BIN=MultiTensor
all: $(BIN)

MultiTensor: MultiTensor.cpp em_update.o tools.o tools_assortative.o cycle_over_realizations.o
	g++ ${CXXFLAGS} ${LIBS2} MultiTensor.cpp ${LIBS} em_update.o tools.o tools_assortative.o cycle_over_realizations.o -o MultiTensor


clean:
	rm -f $(BIN) *.o


em_update.o: em_update.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c em_update.cpp -o em_update.o

tools.o: tools.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c tools.cpp -o tools.o

tools_assortative.o: tools_assortative.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c tools_assortative.cpp -o tools_assortative.o	

cycle_over_realizations.o: cycle_over_realizations.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c cycle_over_realizations.cpp -o cycle_over_realizations.o
