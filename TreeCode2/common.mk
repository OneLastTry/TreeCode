CFLAGS   = -c -Wall -g -I/usr/include/eigen3 -std=c++0x -D_GLIBCXX_PARALLEL
OPTFLAGS = -O3 -msse -msse2 -msse3 -msse4a -ffast-math -march=barcelona -m64 -fopenmp
LDFLAGS  = -g -lboost_system -lboost_unit_test_framework -lboost_chrono -lgslcblas -lsqlite3 -fopenmp 
CC = g++

#Main classes
#SRCS = Node.cpp Particle.cpp Configuration.cpp Tree.cpp
#Boundary conditions
#SRCS += bounds/OpenBoundary.cpp bounds/PeriodicBoundary.cpp
#Potentials
#SRCS += potentials/CoulombForce.cpp potentials/EwaldForce.cpp potentials/InterpolatedEwaldSum.cpp
#Pushers
#SRCS += pushers/LeapfrogPusher.cpp TimeIntegrator.cpp
#Output
#SRCS += output/Output.cpp output/sqlite/DatabaseConnection.cpp

OBJS = $(SRCS:.cpp=.o)

%.o: %.cpp
	$(CC) $(CFLAGS)  $< -o$@