CPP=g++
CFLAGS=-g -Wall
LDFLAGS=-lm
OPTFLAGS=-O3 -ffast-math
MPIFLAGS=-DMPI

NVCC=nvcc
NVCCFLAGS=-DCUDA

PYTHON=python3

serial: main.cpp data.cpp serial.cpp
	$(CPP) $^ -o $@ $(CFLAGS) $(OPTFLAGS)
