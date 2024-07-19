# Makefile for Program TAYLOR

F_COMPILER = gfortran
F_FLAGS    =
F_INCLUDE  = -I/opt/local/include

C_COMPILER = clang++ -O3 -std=c++17
C_FLAGS    = -w 
C_INCLUDE  = -I/opt/local/include 
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas

F_OBJECT = NameListRead.o 
F_SOURCE = NameListRead.f90

C_HEADER = Taylor.h
C_OBJECT = main.o Taylor.o
C_SOURCE = main.cpp Taylor.cpp

EXECUTABLE = ../taylor

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(F_OBJECT) $(EXECUTABLE)

