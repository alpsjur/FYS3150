# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflagsMPI = mpic++ -O3
CPPflags = c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
#LIB = -lboost_system -lboost_filesystem
# Here we define the name of the executable
PROG = main.exe
PROGMPI = mainMPI.exe
TEST_PROG = test_functions.exe

all : target targetMPI tests

target : ${PROG}
targetMPI : ${PROGMPI}
tests : ${TEST_PROG}

${PROG} : main.o metropolis.o ising.o
	${CPPflags} main.o metropolis.o ising.o ${LIB} -o ${PROG}

${PROGMPI} : mainMPI.o metropolis.o ising.o
	${CPPflagsMPI} mainMPI.o metropolis.o ising.o ${LIB} -o ${PROGMPI}

${TEST_PROG} : tests_main.o test_metropolis.o metropolis.o ising.o
	${CPPflags}  tests_main.o test_metropolis.o metropolis.o ising.o ${LIB} -o ${TEST_PROG}


main.o : main.cpp
	${CPPflags} -c main.cpp

mainMPI.o : mainMPI.cpp
	${CPPflagsMPI} -c mainMPI.cpp

metropolis.o : metropolis.cpp
	${CPPflags} -c metropolis.cpp

ising.o : ising.cpp
	${CPPflags} -c ising.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_metropolis.o : test_metropolis.cpp
	${CPPflags} -c test_metropolis.cpp
