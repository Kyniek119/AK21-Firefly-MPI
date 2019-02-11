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
void inicjalizuj_funkcje(int numer_funkcji);
void obliczWartoscFunkcji(int dim, int nf, double* in, double* out);
void zapiszNajlepszyWynik(double* f, double* ffa, double* global_best, double* global_best_param, int ilosc_swietlikow, int wymiar_problemu, double* prev_gen_best);
void obliczKolejnaGeneracje(double* in, double* out, double* val, double alpha, double beta, double gamma, int ilosc_swietlikow, int wymiar_problemu, double prev_gen_best, unsigned int* seed);
void przekopiujObecnaGeneracje(double* new, double* old, int ilosc_swietlikow, int wymiar_problemu);

MPI_Status status;
main(int argc, char **argv){
  unsigned int seed;

  double prev_gen_best = DBL_MAX;
  double global_best = DBL_MAX; //najlepszy wynik globalnie
 
  int i, j, mynum, nprocs;
  int numer_funkcji = 2;

  int ilosc_swietlikow = 500;
  int wymiar_problemu = 300;
  int limit_generacji = 100;
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

  if(mynum==0)startTime = MPI_Wtime();
  for(i=0; i<limit_generacji; i++){
    obliczWartoscFunkcji(
	wymiar_problemu, 
	ilosc_swietlikow, 
	p_ffa, 
	p_f);
    zapiszNajlepszyWynik(
	p_f, 
	p_ffa, 
	&global_best, 
	p_global_best_param, 
	ilosc_swietlikow, 
	wymiar_problemu, 
	&prev_gen_best);
    obliczKolejnaGeneracje(
	p_ffa_prev_gen, 
	p_ffa, 
	p_f, 
	alpha_local, 
	beta_local, 
	gamma_local, 
	ilosc_swietlikow, 
	wymiar_problemu, 
	prev_gen_best, 
	&seed);
    przekopiujObecnaGeneracje(
	p_ffa, 
	p_ffa_prev_gen,
	ilosc_swietlikow, 
	wymiar_problemu); 
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mynum==0){
    endTime = MPI_Wtime();
    printf("Najlepszy wynik: %.4f\n",global_best);
    printf("Czas wykonania obliczeń: %f\n",endTime-startTime);
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

void inicjalizuj_ffa(double* ffa, double* ffa_prev_gen, double* f, int ilosc_swietlikow, int wymiar_problemu){
  int i,j;
  double r;
  double* dane;
  dane = malloc(wymiar_problemu * sizeof(double));
  inicjalizacja_danych(wymiar_problemu, dane);

  for(i=0;i<ilosc_swietlikow;i++){
    for(j=0;j<wymiar_problemu;j++){
      ffa[i*wymiar_problemu+j] = dane[j];
      ffa_prev_gen[i*wymiar_problemu+j] = dane[j];
    }
    f[i] = 1.0;
  }
  free(dane);
}

void obliczWartoscFunkcji(int dim, int nf, double* in, double* out){
  int i, mynum, nprocs;
  double* np;
  MPI_Status status;
  MPI_Request request[nf];
  MPI_Request* p_request = request;
  MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  for(i=mynum; i<nf; i+=nprocs){
    np = &out[i];
    *np = funkcja(dim, &in[i*dim]);
    if(mynum != 0){
      MPI_Isend((void*)&out[i], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, p_request);
    } 
  }

  if(mynum ==0){
    for(i=0; i<nf; i++){
      if(i % nprocs != 0){
        MPI_Recv((void*)&out[i], 1, MPI_DOUBLE, i % nprocs, i, MPI_COMM_WORLD, &status);
      }
    }
  }
  MPI_Bcast((void*)out,nf,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void zapiszNajlepszyWynik(double* f, double* ffa, double* global_best, double* global_best_param, int ilosc_swietlikow, int wymiar_problemu, double* prev_gen_best){
  int mynum, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  *prev_gen_best = DBL_MAX;

  if(mynum == 0){
    int i, id = -1;
    double best = *global_best;
    for( i=0; i<ilosc_swietlikow; i++){
      if(f[i] < best){
        id = i;
        best = f[i];
      }
      if(f[i] < *prev_gen_best){
        *prev_gen_best = f[i];
      }
    }

    if(id != -1){
      //printf("id: %d, with best %.4f \n", id, best);
      *global_best = best;
      memcpy(global_best_param, &ffa[id*wymiar_problemu], sizeof(double) * wymiar_problemu);
    }
  }
  MPI_Bcast(prev_gen_best,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void obliczKolejnaGeneracje(
	double* in, 
	double* out, 
	double* val, 
	double alpha, 
	double beta, 
	double gamma, 
	int ilosc_swietlikow, 
	int wymiar_problemu,
  	double prev_gen_best,
        unsigned int* seed){

  int i, j, k, mynum, nprocs;
  double rnd, beta0, tmp;
  MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  double lokalne_kroki[ilosc_swietlikow * wymiar_problemu];
  double zebranie_krokow[ilosc_swietlikow * wymiar_problemu];
  for(i=0; i<ilosc_swietlikow * wymiar_problemu; i++){
    lokalne_kroki[i] = 0.0;
    zebranie_krokow[i] = 0.0;
  }

  for(i=mynum; i<ilosc_swietlikow; i+=nprocs){
    for(j=0; j<ilosc_swietlikow; j++){
      tmp = 0.0;
      for(k=0;k<wymiar_problemu; k++){
        tmp += (in[i*wymiar_problemu+k] - in[j*wymiar_problemu+k]) * 
		(in[i*wymiar_problemu+k] - in[j*wymiar_problemu+k]);
      }
      tmp = sqrt(tmp);
      if(val[i] > val[j]){
        beta0 = beta * exp(-gamma*pow(tmp,2.0));
        for(k=0; k<wymiar_problemu; k++){
          rnd = ((double)rand_r(seed) / ((double)(RAND_MAX)+(double)(1)));
          double u = alpha * (rnd - 0.5);
          lokalne_kroki[i*wymiar_problemu+k] += beta0 * 
   		(in[j*wymiar_problemu+k] - in[i*wymiar_problemu+k]) + u;
        }
      }
    }
    if(val[i] == prev_gen_best){
        beta0 = beta * exp(-gamma*pow(tmp,2.0));       
        for(k=0; k<wymiar_problemu; k++){
          rnd = ((double)rand_r(seed) / ((double)(RAND_MAX)+(double)(1)));
          double u = alpha * (rnd - 0.5);
          lokalne_kroki[i*wymiar_problemu+k] += u;
        }
    }
  }
  MPI_Allreduce(
	&lokalne_kroki, 
	&zebranie_krokow, 
	ilosc_swietlikow * wymiar_problemu, 
	MPI_DOUBLE,
	MPI_SUM,
	MPI_COMM_WORLD);
  if(mynum==0){
    for(i=0; i<ilosc_swietlikow * wymiar_problemu; i++){
      out[i] += zebranie_krokow[i];
    }
  }
  MPI_Bcast(out,ilosc_swietlikow*wymiar_problemu,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void przekopiujObecnaGeneracje(double* new, double* old, int ilosc_swietlikow, int wymiar_problemu){
  memcpy(old, new, ilosc_swietlikow*wymiar_problemu*sizeof(double));
}
