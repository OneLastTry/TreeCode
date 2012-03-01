LIBDIR = ../../
INCDIR = ../../

CFLAGS   = -c -Wall -g -I/usr/include/eigen3 -std=c++0x -D_GLIBCXX_PARALLEL -I$(INCDIR)
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp
LDFLAGS  = -g -lboost_system -lboost_program_options -lboost_unit_test_framework -lboost_chrono -lgslcblas -lsqlite3 -fopenmp -L$(LIBDIR) 
CC = g++

NAME = test

test: EwaldTest.o
	$(CC) $(LDFLAGS) $< -o $(NAME)
	
clean:
	rm -f EwaldTest.o
	rm -f $(NAME)
	
%.o: %.cpp
	$(CC) $(OPTFLAGS) $(CFLAGS)  $< -o$@ 
