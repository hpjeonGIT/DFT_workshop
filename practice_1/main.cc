#define Nx 2
#define Ny 3
#define Nz 4
#include <iostream>
#include <cmath>
#include <algorithm>
#include "fftw3.h"

double rand_double()
{
  double x = ((double) rand()) / ((double) RAND_MAX);
  return x;
}

int main()
{
  int i, j, k, id, Ntotal;
  Ntotal = Nx*Ny*Nz;

  fftw_complex *ck, *cr, *backup;
  ck = new fftw_complex [Ntotal];
  cr = new fftw_complex [Ntotal];
  backup = new fftw_complex [Ntotal];

  std::cout << "test" << std::endl;
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

  // Backward Transform
  fftw_plan plan;
  plan = fftw_plan_dft_3d(Nx, Ny, Nz, ck, cr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  fftw_destroy_plan(plan);

  // Forward Transform - implement !!!
 

  // Is that all ?

  for (i=0;i<Nx*Ny*Nz;i++) std::cout << ck[i][0] - backup[i][0] << std::endl;

  delete[] backup;
  delete[] ck;
  delete[] cr;
  return 0;
}
