# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = main.exe
PROG2 = body.exe
TEST_PROG = test_functions.exe

all : target target2 tests

target : ${PROG}
target2 : ${PROG2}
tests : ${TEST_PROG}

${PROG} : main.o coordinate.o
	${CPPflags} main.o coordinate.o ${LIB} -o ${PROG}

${PROG2} : body.o coordinate.o
	${CPPflags} body.o coordinate.o ${LIB} -o ${PROG2}

${TEST_PROG} : tests_main.o test_coordinate.o coordinate.o
	${CPPflags} tests_main.o test_coordinate.o coordinate.o ${LIB} -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

coordinate.o : coordinate.cpp
	${CPPflags} -c coordinate.cpp

system.o : system.cpp
	${CPPflags} -c system.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_coordinate.o : test_coordinate.cpp
	${CPPflags} -c test_coordinate.cpp
