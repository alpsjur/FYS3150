# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ #-std=c++0x -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = main.exe
TEST_PROG = test_functions.exe

all : target tests

target : ${PROG}
tests : ${TEST_PROG}

${PROG} : main.o system.o planet.o coordinate.o extractData.o
	${CPPflags} main.o system.o planet.o coordinate.o extractData.o -o ${PROG}

${TEST_PROG} : tests_main.o test_coordinate.o coordinate.o
	${CPPflags} tests_main.o test_coordinate.o coordinate.o -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

coordinate.o : coordinate.cpp
	${CPPflags} -c coordinate.cpp

planet.o : planet.cpp
	${CPPflags} -c planet.cpp

system.o : system.cpp
	${CPPflags} -c system.cpp

extractData.o : extractData.cpp
	${CPPflags} -c extractData.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_coordinate.o : test_coordinate.cpp
	${CPPflags} -c test_coordinate.cpp
