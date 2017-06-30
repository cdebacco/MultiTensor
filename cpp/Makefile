#CXXFLAGS=-Wall -O0 -g -Wno-deprecated -lgomp -fopenmp
CXXFLAGS=-Wall -std=c++0x -O3 -g -Wno-deprecated
#LIBS=-L/usr/local/Cellar/boost/1.58.0/lib/ -lboost_program_options
LIBS=-L/usr/local/Cellar/boost/1.64.0_1/lib -lboost_program_options
LIBS2=-I /usr/include
BIN=MultiTensor_undirected
all: $(BIN)

MultiTensor: MultiTensor.cpp em_update.o tools.o tools_assortative.o cycle_over_realizations.o
	g++ ${CXXFLAGS} ${LIBS2} MultiTensor.cpp ${LIBS} em_update.o tools.o tools_assortative.o cycle_over_realizations.o -o MultiTensor
MultiTensor_undirected: MultiTensor_undirected.cpp em_update_undirected.o tools.o tools_assortative.o cycle_over_realizations_undirected.o
	g++ ${CXXFLAGS} ${LIBS2} MultiTensor_undirected.cpp ${LIBS} em_update_undirected.o tools.o tools_assortative.o cycle_over_realizations_undirected.o -o MultiTensor_undirected


clean:
	rm -f $(BIN) *.o


em_update.o: em_update.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c em_update.cpp -o em_update.o

em_update_undirected.o: em_update_undirected.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c em_update_undirected.cpp -o em_update_undirected.o

tools.o: tools.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c tools.cpp -o tools.o

tools_assortative.o: tools_assortative.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c tools_assortative.cpp -o tools_assortative.o	

cycle_over_realizations.o: cycle_over_realizations.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c cycle_over_realizations.cpp -o cycle_over_realizations.o
cycle_over_realizations_undirected.o: cycle_over_realizations_undirected.cpp mlg.hpp 
	g++ ${CXXFLAGS} -c cycle_over_realizations_undirected.cpp -o cycle_over_realizations_undirected.o

em_uni_dis_r: em_uni_dis_r.cpp
	g++ ${CXXFLAGS} ${LIBS2} em_uni_dis_r.cpp ${LIBS} -o em_uni_dis_r

em_uni_as: em_uni_as.cpp
	g++ ${CXXFLAGS} ${LIBS2} em_uni_as.cpp ${LIBS} -o em_uni_as