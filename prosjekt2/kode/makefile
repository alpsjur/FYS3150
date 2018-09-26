# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG=jacobis_method.exe
${PROG} :	   	jacobis_method.o
			${CPPflags} jacobis_method.o ${LIB} -o ${PROG}

jacobis_method.o :		jacobis_method.cpp
		        	${CPPflags} -c jacobis_method.cpp
