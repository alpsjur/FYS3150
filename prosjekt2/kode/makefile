# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = jacobis_method.exe
TEST_PROG = test_jacobi.exe

all : target tests

target : ${PROG}
tests : ${TEST_PROG}

${PROG} :	   	jacobis_method.o
			${CPPflags} jacobis_method.o ${LIB} -o ${PROG}

${TEST_PROG} :			test_jacobi.o jacobis_method.o
			${CPPflags} test_jacobi.o jacobis_method.o ${LIB} -o ${TEST_PROG}

jacobis_method.o :		jacobis_method.cpp
		        	${CPPflags} -c jacobis_method.cpp


test_jacobi.o : 			test_jacobi.cpp
							${CPPflags} -c $<
