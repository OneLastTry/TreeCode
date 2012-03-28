INCDIR = ../../

CFLAGS   = -Wall -g -I/usr/include/eigen3 -D_GLIBCXX_PARALLEL -I$(INCDIR)
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp
LDFLAGS  = -g -lboost_system -lboost_program_options -lboost_unit_test_framework -fopenmp 
CC = gfilt

NAME = ewald 

all: 
	$(CC) $(CFLAGS) $(LDFLAGS) $(OPTFLAGS) main.cpp -o $(NAME)
	
clean:
	rm -f $(NAME)
	
