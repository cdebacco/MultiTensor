#CXXFLAGS=-Wall -O0 -g -Wno-deprecated -lgomp -fopenmp
CXXFLAGS=-Wall -std=c++0x -O0 -g -Wno-deprecated
LIBS=-L/usr/local/Cellar/boost/1.58.0/lib/ -lboost_program_options
LIBS2=-I /usr/include
BIN=em_dis
all: $(BIN)
em: em.cpp 
	g++ ${CXXFLAGS} ${LIBS2} em.cpp ${LIBS} -o em

em_dis: em_dis.cpp 
	g++ ${CXXFLAGS} ${LIBS2} em_dis.cpp ${LIBS} -o em_dis

em_uni_as: em_uni_as.cpp 
	g++ ${CXXFLAGS} ${LIBS2} em_uni_as.cpp ${LIBS} -o em_uni_as

em_uni_dis: em_uni_dis.cpp 
	g++ ${CXXFLAGS} ${LIBS2} em_uni_dis.cpp ${LIBS} -o em_uni_dis


clean:
	rm -f $(BIN) *.o



