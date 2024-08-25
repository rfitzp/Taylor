# Makefile for Program TAYLOR

C_COMPILER = clang++ 
C_FLAGS    = -w -O3 -std=c++17 
C_INCLUDE  = -I/opt/local/include 
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgsl -lgslcblas

C_HEADER = Taylor.h
C_SOURCE = main.cpp Taylor.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../taylor

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

