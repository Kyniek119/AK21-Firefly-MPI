#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funkcje.h"
#include <float.h>
#include "mpi.h"
#define MAX_D 1000
#define MAX_FFA 1000

double ffa[MAX_FFA][MAX_D]; //swietliki
double ffa_prev_gen[MAX_FFA][MAX_D]; //zmienna do przechowywania poprzedniej generacji.
double f[MAX_FFA]; //wartosc funkcji


double global_best = DBL_MAX; //najlepszy wynik globalnie
double global_best_param[MAX_D];

typedef double (*FunctionalCallback)(int dimension, double sol[MAX_D]);
FunctionalCallback funkcja = &funkcja1;

typedef void (*InitializationCallback)(int dimension, double sol[MAX_D]);
InitializationCallback inicjalizacja_danych = &init1;

void inicjalizuj_ffa(int ilosc_swietlikow, int wymiar_problemu);
void inicjalizuj_funkcje(int numer_funkcji);
void obliczWartosciFunkcji(int dim, double* in, double* out);
void zapiszNajlepszyWynik(double f[MAX_FFA], double ffa[MAX_FFA][MAX_D], double* global_best, double* global_best_param[MAX_D], int ilosc_swietlikow);

MPI_Status status;
main(int argc, char **argv)
{
float err, sum, w, x;
int i, j, N, info, step, mynum, nprocs, source, dest = 0;
int type = 2, nbytes = 0, EUI_SUCCEED = 0;
int numer_funkcji = 1;

int ilosc_swietlikow = 10;
int wymiar_problemu = 20;
int limit_generacji = 2;
double alpha_local = 0.5;
double beta_local = 1;
double gamma_local = 0.01;

MPI_Request request[ilosc_swietlikow];
MPI_Status status[ilosc_swietlikow];
int requestNumber = 0;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

inicjalizuj_funkcje(numer_funkcji);
inicjalizuj_ffa(ilosc_swietlikow, wymiar_problemu);

step = (int)(ilosc_swietlikow / nprocs);

for (i = 0; i < limit_generacji; i++){
  requestNumber = 0;
  for(j = mynum; j < ilosc_swietlikow; j += nprocs){
    obliczWartosciFunkcji(wymiar_problemu, ffa[j], &f[j]);
    if(mynum != 0){
      info = MPI_Isend((void*)&f[j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request[requestNumber++]);
      printf("Wysłano wartosc %.4f z watku %d\n", f[j], mynum);
    }
  }

  if(mynum == 0){
    for(j = 0; j < ilosc_swietlikow; j++){
      if(j % nprocs != 0){
        info = MPI_Irecv((void*)&f[j], 1, MPI_DOUBLE, j % nprocs, 0, MPI_COMM_WORLD, &request[requestNumber++]);
      }
    }
  }

  for(j=0; j < requestNumber; j++){
    MPI_Wait(&request[j], &status[j]);
  }
  
  if(mynum == 0){
    zapiszNajlepszyWynik(f, ffa, &global_best, &global_best_param, ilosc_swietlikow);
  }
 
  //obliczKolejnaGeneracje(in, out, alpha_local, beta_local, gamma_local);
  
}
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

void inicjalizuj_ffa(int ilosc_swietlikow, int wymiar_problemu){
  int i,j;
  double r;
  double dane[MAX_D];
  inicjalizacja_danych(wymiar_problemu, dane);

  for(i=0;i<ilosc_swietlikow;i++){
    for(j=0;j<wymiar_problemu;j++){
      ffa[i][j] = dane[j];
      ffa_prev_gen[i][j] = dane[j];
    }
    f[i] = 1.0;
  }
}

void obliczWartosciFunkcji(int dim, double* in, double* out){
  double res = funkcja(dim, in);
  *out = res;
}

void zapiszNajlepszyWynik(double f[MAX_FFA], double ffa[MAX_FFA][MAX_D], double* global_best, double* global_best_param[MAX_D], int ilosc_swietlikow){
  int i, id = -1;
  double best = *global_best;
  for( i=0; i<ilosc_swietlikow; i++){
    if(f[i] < best){
      id = i;
      best = f[i];
    }
  }

  if(id != -1){
    printf("id: %d, with best %.4f \n", id, best);
    *global_best = best;
    memcpy(global_best_param, ffa[id], sizeof(double) * MAX_D);
  }
}
