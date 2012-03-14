LIBDIR = ../../
INCDIR = ../../

CFLAGS   = -Wall -g -I/usr/include/eigen3 -std=c++0x -D_GLIBCXX_PARALLEL -I$(INCDIR)
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp
LDFLAGS  = -g -lboost_program_options -fopenmp -L$(LIBDIR)
CC = g++

NAME = resistivity-calculator

all:
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTFLAGS) main.cpp -o $(NAME)
	
clean:
	rm -f $(NAME)
 
