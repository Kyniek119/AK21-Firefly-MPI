CC=mpicc.openmpi
GCC=gcc
CFLAGS=-O2 -lm
DEPS = Firefly_MPI.h Funkcje.h

%.o: %.c $(DEPS)
	$(GCC) -o $@ -lm $<

FireflyMPI: Firefly_MPI.o Funkcje.o
	$(CC) -o Firefly_MPI.out Funkcje.o Firefly_MPI.o $(CFLAGS) 
