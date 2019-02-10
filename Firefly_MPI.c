#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Funkcje.h"
#include <float.h>
#include "mpi.h"



typedef double (*FunctionalCallback)(int dimension, double* sol);
FunctionalCallback funkcja = &funkcja1;

typedef void (*InitializationCallback)(int dimension, double* sol);
InitializationCallback inicjalizacja_danych = &init1;

void inicjalizuj_ffa(double* ffa, double* ffa_prev_gen, double* f, int ilosc_swietlikow, int wymiar_problemu);
void cleanUp(double* ffa, double* ffa_prev_gen, double* f, double* global_best_param);
void inicjalizuj_funkcje(int numer_funkcji);
void obliczWartoscFunkcji(int dim, int nf, double* in, double* out);
void zapiszNajlepszyWynik(double* f, double* ffa, double* global_best, double* global_best_param, int ilosc_swietlikow, int wymiar_problemu);
void obliczKolejnaGeneracje(double* in, double* out, double val, double alpha, double beta, double gamma, int ilosc_swietlikow, int wymiar_problemu, int mynum, int nprocs);

MPI_Status status;
main(int argc, char **argv){
  unsigned int seed;

  double global_best = DBL_MAX; //najlepszy wynik globalnie
 
  int i, j, mynum, nprocs;
  int numer_funkcji = 1;

  int ilosc_swietlikow = 10;
  int wymiar_problemu = 20;
  int limit_generacji = 2;
  double alpha_local = 0.5;
  double beta_local = 1;
  double gamma_local = 0.01;

  double ffa[ilosc_swietlikow*wymiar_problemu]; //swietliki
  double ffa_prev_gen[ilosc_swietlikow*wymiar_problemu]; //zmienna do przechowywania poprzedniej generacji.
  double f[ilosc_swietlikow]; //wartosc funkcji
  double global_best_param[wymiar_problemu];

  double* p_ffa = ffa; //swietliki
  double* p_ffa_prev_gen = ffa_prev_gen; //zmienna do przechowywania poprzedniej generacji.
  double* p_f = f; //wartosc funkcji
  double* p_global_best_param = global_best_param;

  double startTime, endTime;

  MPI_Request request[ilosc_swietlikow];
  MPI_Status status[ilosc_swietlikow];
  int requestNumber = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  inicjalizuj_funkcje(numer_funkcji);
  inicjalizuj_ffa(p_ffa, p_ffa_prev_gen, p_f, ilosc_swietlikow, wymiar_problemu);
  seed = time(NULL) ^ getpid();

  startTime = MPI_Wtime();
  for(i=0; i<limit_generacji; i++){
    obliczWartoscFunkcji(wymiar_problemu, ilosc_swietlikow, p_ffa, p_f);
    //sendCalculationResults();
    //zapiszNajlepszyWynik();
    //obliczKolejnaGeneracje();
    //sendNewlyCalculatedGeneration();    
  }
  endTime = MPI_Wtime();
  printf("Czas wykonania obliczeń: %f\n",endTime-startTime);

  cleanUp(ffa, ffa_prev_gen, f, global_best_param);
  MPI_Finalize();

}

void inicjalizuj_funkcje(int numer_funkcji){
  switch(numer_funkcji){
	case 1: funkcja = &funkcja1;
		inicjalizacja_danych = &init1;
		break;
        case 2: funkcja = &funkcja2;
		inicjalizacja_danych = &init2;
		break;
        case 3: funkcja = &funkcja3;
		inicjalizacja_danych = &init3;
		break;
	default: funkcja = &funkcja1;
                 inicjalizacja_danych = &init1;
		break;
  }
} 

void inicjalizuj_ffa(double* ffa, double* ffa_prev_gen, double* f, int ilosc_swietlikow, int wymiar_problemu){
  int i,j;
  double r;
  double* dane;
  dane = malloc(wymiar_problemu * sizeof(double));
  inicjalizacja_danych(wymiar_problemu, dane);

  for(i=0;i<ilosc_swietlikow;i++){
    for(j=0;j<wymiar_problemu;j++){
      *(ffa+i*wymiar_problemu+j) = *(dane+j);
      *(ffa_prev_gen+i*wymiar_problemu+j) = *(dane+j);
    }
    *(f+i) = 1.0;
    printf("Podstawiono %f pod wartosc f[%d], %p\n",*(f+i), i, f+i);
  }

  free(dane);
}

void obliczWartoscFunkcji(int dim, int nf, double* in, double* out){
  int i, mynum, nprocs;
  MPI_Request request[nf];
  MPI_Request* p_request = request;
  MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  for(i=mynum; i<nf; i+=nprocs){
    *(out+i) = funkcja(dim, in+i*dim);
    if(mynum != 0){
      //MPI_Isend(out+i, 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, p_request);
    } 
  }

  if(mynum ==0){
    for(i=0; i<nf; i++){
      if(i % mynum != 0){
        //MPI_Irecv(out+i, 1, MPI_DOUBLE, i % mynum, i, MPI_COMM_WORLD, p_request+i);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
}

void zapiszNajlepszyWynik(double* f, double* ffa, double* global_best, double* global_best_param, int ilosc_swietlikow, int wymiar_problemu){
  int i, id = -1;
  double best = *global_best;
  for( i=0; i<ilosc_swietlikow; i++){
    if(*(f+i) < best){
      id = i;
      best = *(f+i);
    }
  }

  if(id != -1){
    printf("id: %d, with best %.4f \n", id, best);
    *global_best = best;
    memcpy(global_best_param, ffa+id*wymiar_problemu, sizeof(double) * wymiar_problemu);
  }
}

void cleanUp(double* ffa, double* ffa_prev_gen, double* f, double* global_best_param){
  
}
