# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3
# Here we define the library functions we nee
LIB = #-larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = rossby.exe
TEST_PROG = test_jacobi.exe

all : target

target : ${PROG}
tests : ${TEST_PROG}

${PROG} : main.o rossby.o
	${CPPflags} main.o rossby.o ${LIB} -o ${PROG}

${TEST_PROG} : tests_main.o test_jacobi.o jacobis_method.o
	${CPPflags} tests_main.o test_jacobi.o jacobis_method.o ${LIB} -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

jacobis_method.o : rossby.cpp
	${CPPflags} -c rossby.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_jacobi.o : test_jacobi.cpp
	${CPPflags} -c test_jacobi.cpp
