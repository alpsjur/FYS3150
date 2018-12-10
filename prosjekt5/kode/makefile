# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROGPERIODIC = rossby_periodic.exe
PROGBOUNDED = rossby_bounded.exe
TEST_PROG = test_jacobi.exe

all : king kong

king : ${PROGPERIODIC}
kong : ${PROGBOUNDED}
tests : ${TEST_PROG}

${PROGPERIODIC} : rossby_periodic.o
	${CPPflags} rossby_periodic.o ${LIB} -o ${PROGPERIODIC}

${PROGBOUNDED} : rossby_bounded.o
	${CPPflags} rossby_bounded.o -o ${PROGBOUNDED}

${TEST_PROG} : tests_main.o
	${CPPflags} tests_main.o ${LIB} -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

rossby_periodic.o : rossby_periodic.cpp
	${CPPflags} -c rossby_periodic.cpp

rossby_bounded.o : rossby_bounded.cpp
	${CPPflags} -c rossby_bounded.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp
