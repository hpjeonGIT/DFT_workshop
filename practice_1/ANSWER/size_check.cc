#define Nx 128
#define Ny 128
#define Nz 128
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include "fftw3.h"

double rand_double()
{
  double x = ((double) rand()) / ((double) RAND_MAX);
  return x;
}

int main()
{
 time_t t0, t1;
  t0 = time(NULL);
  int i, j, k, id, Ntotal;
  Ntotal = Nx*Ny*Nz;

  fftw_complex *ck, *cr, *backup;
  ck = new fftw_complex [Ntotal];
  cr = new fftw_complex [Ntotal];
  backup = new fftw_complex [Ntotal];

  //std::cout << "test" << std::endl;
  for(i=0;i<Nx;i++) { 
    for (j=0;j<Ny;j++) {
      for (k=0;k<Nz;k++) {
	id = i + Nx*(j+Ny*k);
	ck[id][0] = rand_double();
	ck[id][1] = rand_double();
	backup[id][0] = ck[id][0];
	backup[id][1] = ck[id][1];
      }
    }
  }
  for (int n = 0;n<20;n++) {
    fftw_plan plan;
    plan = fftw_plan_dft_3d(Nx, Ny, Nz, ck, cr, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan); 
    
    fftw_destroy_plan(plan);
    plan = fftw_plan_dft_3d(Nx, Ny, Nz, cr, ck, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan); 
    fftw_destroy_plan(plan);
 
    for(i=0;i<Nx*Ny*Nz;i++) {
      ck[i][0] = ck[i][0]/Ntotal;
      ck[i][1] = ck[i][1]/Ntotal;
    }
  }
  //for (i=0;i<Nx*Ny*Nz;i++) std::cout << ck[i][0] - backup[i][0] << std::endl;
  
  delete[] backup;
  delete[] ck;
  delete[] cr;
  t1 = time(NULL);
  std::cout << t0 << " " << t1 << " "  << t1 - t0 << std::endl;
  return 0;
}