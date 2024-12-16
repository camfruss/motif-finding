CPP=g++ -std=c++20 
CFLAGS=-lm -g -Wall
OPTFLAGS=-O3 -ffast-math
MPIFLAGS=-DMPI

NVCC=nvcc
NVCCFLAGS=-DCUDA

PYTHON=python3

SOURCES=main.cpp data.cpp serial.cpp utility.cpp
OBJECTS=$(SOURCES:.cpp=.o)
DEPS=data.hpp serial.hpp utility.hpp

TARGETS=serial

all: $(TARGETS)

serial: $(SOURCES) 
	$(CPP) $^ -o $@ $(CFLAGS) $(OPTFLAGS)

clean:
	rm -f $(OBJECTS) $(TARGETS)
