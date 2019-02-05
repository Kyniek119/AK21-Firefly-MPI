#include <math.h>
#include "Funkcje.h"

//funkcja "init1" inicjalizuje tablice wartosci dla funkcji "funkcja1", jest ona zapisywana w dostarczonym buforze "sol".
void init1(int wymiar_problemu, double* sol){
  int i;
  for(i=0;i<wymiar_problemu;i++){
    sol[i] = 3.0;
  }
}

//implementacja "1. Quadratic function" z funkcji testowych do zadań optymalizacji.
double funkcja1(int wymiar_problemu, double* sol){
  double sum = 0.0;
  int i;
  for(i=2;i<wymiar_problemu;i++){
    sum += 	100*(pow(sol[i],2.0)+pow(sol[i-1],2.0)) 
		+pow(sol[i-2],2.0); 
  }
  return sum;
}

//funkcja "init2" inicjalizuje tablice wartosci dla funkcji "funkcja2", jest ona zapisywana w dostarczonym buforze "sol".
void init2(int wymiar_problemu, double* sol){
  int i;
  for(i=0;i<wymiar_problemu;i++){
    if(i%2 == 1) sol[i] = -1.0;
    else sol[i] = -3.0;
  }
}

//implementacja "2. Woods function" z funkcji testowych do zadań optymalizacji.
double funkcja2(int wymiar_problemu, double* sol){
  double sum = 0.0;
  int i;
  int suma_do = (int)(wymiar_problemu/4);
  for(i=1;i<=suma_do;i++){
    sum+= 	100*pow((sol[4*i-2]-pow(sol[4*i-3],2.0)),2.0)
		+pow((1-sol[4*i-3]),2.0)
		+90*pow((sol[4*i]-pow(sol[4*i-1],2.0)),2.0)
		+pow((1-sol[4*i-1]),2.0)
		+10*pow((sol[4*i-2]+sol[4*i]-2),2.0)
		+0.1*pow((sol[4*i-2]-sol[4*i]),2.0); 
  }
  return sum;
}

//funkcja "init3" inicjalizuje tablice wartosci dla funkcji "funkcja3", jest ona zapisywana w dostarczonym buforze "sol".
void init3(int wymiar_problemu, double* sol){
  int i;
  for(i=0;i<wymiar_problemu;i++){
    if(i%4 == 1) sol[i] = -1.0;
    else if (i%4 == 2) sol[i] = 0.0;
    else if (i%4 == 3) sol[i] = 1.0;
    else sol[i] = 3.0;
  }
}

//implementacja "3. Powell singular function" z funkcji testowych do zadań optymalizacji.
double funkcja3(int wymiar_problemu, double* sol){
  double sum = 0.0;
  int i;
  int suma_do = (int)(wymiar_problemu/4);
  for(i=1;i<=suma_do;i++){
    sum+= 	pow((sol[4*i-3]+10*sol[4*i-2]), 2.0)
		+5*pow((sol[4*i-1]-sol[4*i]),2.0)
		+pow(sol[4*i-2]-2*sol[4*i-1],4.0)
		+10*pow(sol[4*i-3]-sol[4*i],4.0); 
  }
  return sum;
}
