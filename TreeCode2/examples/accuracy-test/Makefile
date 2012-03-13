INCDIR = ../../

CFLAGS   = -Wall -g -I/usr/include/eigen3 -std=c++0x -I$(INCDIR)
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp -D_GLIBCXX_PARALLEL 
LDFLAGS  = -fopenmp 
CC = g++

NAME = accuracy-test

all:
	$(CC) $(LDFLAGS) $(CFLAGS) $(OPTFLAGS) main.cpp -o $(NAME)
	
clean:
	rm -f $(NAME)
	
