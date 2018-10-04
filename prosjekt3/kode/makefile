# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = main.exe
TEST_PROG = test_functions.exe

all : target
#tests

target : ${PROG}
#tests : ${TEST_PROG}

${PROG} : main.o coordinate.o
	${CPPflags} main.o coordinate.o ${LIB} -o ${PROG}

#${TEST_PROG} : tests_main.o test_jacobi.o jacobis_method.o
	#${CPPflags} tests_main.o test_jacobi.o jacobis_method.o ${LIB} -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

coordinate.o : coordinate.cpp
	${CPPflags} -c coordinate.cpp
