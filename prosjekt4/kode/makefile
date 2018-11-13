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
TEST_PROG = test_functions.exe

all : target tests

target : ${PROG}
tests : ${TEST_PROG}

${PROG} : main.o metropolis.o ising.o
	${CPPflagsMPI} main.o metropolis.o ising.o ${LIB} -o ${PROG}

${TEST_PROG} : tests_main.o test_metropolis.o metropolis.o ising.o
	${CPPflags}  tests_main.o test_metropolis.o metropolis.o ising.o ${LIB} -o ${TEST_PROG}


main.o : main.cpp
	${CPPflagsMPI} -c main.cpp

metropolis.o : metropolis.cpp
	${CPPflagsMPI} -c metropolis.cpp

ising.o : ising.cpp
	${CPPflagsMPI} -c ising.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_metropolis.o : test_metropolis.cpp
	${CPPflags} -c test_metropolis.cpp
