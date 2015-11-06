#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"



extern "C" void zheev_(char *job, char *uplo, int *n,  dcomplex *a, int *lda, 
                       double *w, dcomplex *work,  int *lwork, double *rwork, 
                       int *info);

//////////////////////////////////////////////////////////////////////
void Hermitian_solv(int pw[][3], fftw_ **Hht, fftw_ **Hloc, fftw_ *Vxc,
		    double *Hkin, fftw_ *ck)
{

  int i, j, n, ix, iy, iz, jx, jy, jz, kx, ky, kz, id, id2;

  double *w=NULL, *rwork=NULL;
  dcomplex *a=NULL, *work=NULL;
  int info=0, lwork;
  a = new dcomplex [Npw*Npw];
  w = new double [Npw];
  rwork = new double [3*Npw-2];
  int nspan = Npw;
  for (i=0;i<Npw;i++) {
    ix = pw[i][0]; iy = pw[i][1]; iz = pw[i][2];
    for (j=0;j<Npw;j++) {
      jx = pw[j][0]; jy = pw[j][1]; jz = pw[j][2];
      kx = ix - jx;kx = (kx+Ngr)%Ngr;
      ky = iy - jy;ky = (ky+Ngr)%Ngr;
      kz = iz - jz;kz = (kz+Ngr)%Ngr;
      id = kz + Ngr*(ky + Ngr*kx);
      
      a[i+Npw*j].real = // Fill out
      a[i+Npw*j].imag = // Fill out
      if (i == j)  a[i+Npw*j].real += Hkin[i];
    }
  }
  work = new dcomplex [1];  
  lwork = -1;
  // query
  zheev_( "V", "U", &nspan, a, &nspan, w, work, &lwork, rwork, &info );
  lwork = (int) work[0].real;
  delete[] work;
  work = new dcomplex [lwork];

  // actual solving
  zheev_( "V", "U", &nspan, a, &nspan, w, work, &lwork, rwork, &info );
  // w - eigen value   If INFO = 0, the eigenvalues in ascending order.


  std::cout << "Min energy is " << w[0] << " "  << w[1] << " " << w[2] 
            << " " << w[3] <<  " " << w[4] <<std::endl;

  for (i=0;i<Npw;i++) {
    ix = pw[i][0]; iy = pw[i][1]; iz = pw[i][2];
    ix = (ix+Ngr)%Ngr; iy = (iy+Ngr)%Ngr; iz = (iz+Ngr)%Ngr;
    id = iz + Ngr*(iy + Ngr*ix);
    id2 = i + 0*Npw;
    ck[id][0] = a[id2].real;
    ck[id][1] = a[id2].imag;
  }
  delete[] w; delete[] rwork; delete[] a; delete[] work;
}
