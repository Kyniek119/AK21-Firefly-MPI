CC=mpicc.openmpi
CFLAGS=-O2 -lm
DEPS = Firefly_MPI.h Funkcje.h

%.o: %.c $(DEPS)
	$(CC) -o $@ -lm $<

FireflyMPI: Firefly_MPI.o Funkcje.o
	$(CC) -o Firefly_MPI.out Funkcje.o Firefly_MPI.o $(CFLAGS) 
