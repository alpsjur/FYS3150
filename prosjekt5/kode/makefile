# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROGPERIODIC = rossby_periodic.exe
PROGBOUNDED = rossby_bounded.exe
PROGPERIODIC2D = rossby_periodic_2d.exe
PROGBOUNDED2D = rossby_bounded_2d.exe
TEST_PROG = test_jacobi.exe

all : king kong ali baba

king : ${PROGPERIODIC}
kong : ${PROGBOUNDED}
ali : ${PROGPERIODIC2D}
baba : ${PROGBOUNDED2D}
tests : ${TEST_PROG}

${PROGPERIODIC} : rossby_periodic.o
	${CPPflags} rossby_periodic.o ${LIB} -o ${PROGPERIODIC}

${PROGBOUNDED} : rossby_bounded.o
	${CPPflags} rossby_bounded.o -o ${PROGBOUNDED}

${PROGPERIODIC2D} : rossby_periodic_2d.o
	${CPPflags} rossby_periodic_2d.o ${LIB} -o ${PROGPERIODIC2D}

${PROGBOUNDED2D} : rossby_bounded_2d.o
	${CPPflags} rossby_bounded_2d.o ${LIB} -o ${PROGBOUNDED2D}

${TEST_PROG} : tests_main.o
	${CPPflags} tests_main.o ${LIB} -o ${TEST_PROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

rossby_periodic.o : rossby_periodic.cpp
	${CPPflags} -c rossby_periodic.cpp

rossby_bounded.o : rossby_bounded.cpp
	${CPPflags} -c rossby_bounded.cpp

rossby_periodic_2d.o : rossby_periodic_2d.cpp
	${CPPflags} -c rossby_periodic_2d.cpp

rossby_bounded_2d.o : rossby_bounded_2d.cpp
	${CPPflags} -c rossby_bounded_2d.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp
