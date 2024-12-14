CPP=g++
CFLAGS=
OPTFLAGS=-O3 -ffast-math
MPIFLAGS=-DMPI

NVCC=nvcc
NVCCFLAGS=-DCUDA

PYTHON=python3

serial: main.cpp data.cpp
	$(CPP) $^ $@ $(CFLAGS) $(OPTFLAGS)
