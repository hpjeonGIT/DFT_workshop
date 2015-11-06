#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"

//
// Random number generator
double rand_double() // random number generator
{
  double x = (double)(rand())/(double)(RAND_MAX) - 0.5;
  return x;
}

//
// Initialization of reciprocal wave vectors using random numbers
void ck_initialize (int pw[][3], fftw_ *ck, fftw_ *nr, fftw_ *nk)
{
  int i, kx, ky, kz, id;
  double x;
  int Ngr3 = Ngr*Ngr*Ngr;
  
  for (i=0;i<Ngr3;i++) ck[i][0] = ck[i][1] = 0.0;
  for (i=0;i<Ngr3;i++) nr[i][0] = nr[i][1] = nk[i][0] = nk[i][1] = 0.0;
  
  //
  // Initialization and normalization
  x = 0.0;
  for (i=0;i<Npw;i++) {
    kx = (pw[i][0] + Ngr)%Ngr;
    ky = (pw[i][1] + Ngr)%Ngr;
    kz = (pw[i][2] + Ngr)%Ngr;
    id = kz + Ngr*(ky + Ngr*kx);
    ck[id][0] =  rand_double();
    ck[id][1] =  rand_double();
    x += ck[id][0]*ck[id][0] + ck[id][1]*ck[id][1];
  }
  x = sqrt(x);
  for (i=0;i<Ngr3;i++) {
    ck[i][0] = ck[i][0]/x;
    ck[i][1] = ck[i][1]/x;
  }
}

//
// |phi> => <r|phi>
void ck_fft_cr(fftw_ *ck, fftw_ *cr, double l_box)
{
  int i;
  int Ngr3 = Ngr*Ngr*Ngr;
  fftw_plan plan;
  plan = fftw_plan_dft_3d(Ngr, Ngr, Ngr, ck, cr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan); fftw_destroy_plan(plan); 
  
  double omega = pow(l_box, 3.);

  for (i=0;i<Ngr3;i++) {
    cr[i][0] = cr[i][0]/sqrt(omega);
    cr[i][1] = cr[i][1]/sqrt(omega);
  } 
}
				
//
//
void nr_fft_nk(fftw_ *cr, fftw_ *nr, fftw_ *nk, double alpha)
{
  int i;
  int Ngr3 = Ngr*Ngr*Ngr;
  double rho[Ngr3];

  for(i=0;i<Ngr3;i++) {      
    rho[i] = cr[i][0]*cr[i][0] + cr[i][1]*cr[i][1];
  }

  for (i=0;i<Ngr3;i++) {
    nr[i][0] = (1. - alpha)*nr[i][0] + alpha*rho[i];
  }
  fftw_plan plan;
  plan = fftw_plan_dft_3d(Ngr, Ngr, Ngr, nr, nk, FFTW_FORWARD, FFTW_ESTIMATE); 
  fftw_execute(plan); fftw_destroy_plan(plan); 
  for (i=0;i<Ngr3;i++) {
    nk[i][0] = nk[i][0]/(double) Ngr3;
    nk[i][1] = nk[i][1]/(double) Ngr3;
  }
}
