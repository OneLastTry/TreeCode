LIBDIR = ../../
INCDIR = ../../

CFLAGS   = -c -Wall -g -I/usr/include/eigen3 -std=c++0x -D_GLIBCXX_PARALLEL -I$(INCDIR)
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp
LDFLAGS  = -g -lboost_system -lboost_program_options -lboost_unit_test_framework -lboost_chrono -lgslcblas -lsqlite3 -fopenmp -L$(LIBDIR) -ltreecode 
CC = g++

NAME = resistivity-calculator

all: main.o
	$(CC) $(LDFLAGS) main.o -o $(NAME)
	
clean:
	rm -f main.o
	rm -f $(NAME)
	
%.o: %.cpp
	$(CC) $(OPTFLAGS) $(CFLAGS)  $< -o$@ 
