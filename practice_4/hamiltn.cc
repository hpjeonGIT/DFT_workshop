#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"

void V_HT(int pw[][3], fftw_ *nk, fftw_ **Hht, double &Eht, double l_box)
{
  int i, j, ix, iy, iz, jx, jy, jz, kx, ky, kz, id;
  double G2;
  double omega = pow(l_box, 3.);

  Eht = 0.0;
 
  for (i=0;i<Npw;i++) {
    ix = pw[i][0]; iy = pw[i][1]; iz = pw[i][2];
    for (j=0;j<Npw;j++) {
      jx = pw[j][0]; jy = pw[j][1]; jz = pw[j][2];
      kx = ix - jx; ky = iy - jy; kz = iz - jz;
      if (i != j) {
	G2 = FOURPISQ * (double) (kx*kx + ky*ky + kz*kz) / (l_box*l_box);
	kx = (kx + Ngr)%Ngr; ky = (ky + Ngr)%Ngr; kz = (kz + Ngr)%Ngr;
	id = kz + Ngr*(ky + Ngr*kx);
	Hht[i][j][0] = // Fill out
	Hht[i][j][1] = // Fill out
	Eht += // Fill out
      }
      else {
	Hht[i][j][0] = Hht[i][j][1] = 0.0;
      }
    }
  }
}

//
//
void V_loc(int pw[][3], fftw_ *nk, fftw_ **Hloc, double &Eloc, double l_box)
{
  int i, j, ix, iy, iz, jx, jy, jz, kx, ky, kz, id;
  double G2;
  double omega = pow(l_box, 3.);
  double rs, c1, c2, A;
  rs = 0.2;
  c1 = -4.0663326;
  c2 = 0.6778322;
  
  Eloc = 0.0;

  for (i=0;i<Npw;i++) {
    ix = pw[i][0]; iy = pw[i][1]; iz = pw[i][2];
    for (j=0;j<Npw;j++) {
      jx = pw[j][0]; jy = pw[j][1]; jz = pw[j][2];
      kx = ix - jx; ky = iy - jy; kz = iz - jz;
      if (i != j) {
	G2 = FOURPISQ * (double) (kx*kx + ky*ky + kz*kz) / (l_box*l_box);
	kx = (kx + Ngr)%Ngr; ky = (ky + Ngr)%Ngr; kz = (kz + Ngr)%Ngr;
	id = kz + Ngr*(ky + Ngr*kx);
	A = exp(-G2*rs*rs*0.5);
	Hloc[i][j][0] = // Fill out
	Hloc[i][j][1] = 0.0;
	Eloc += // Fill out 
      }
      else {
	Hloc[i][j][0] = Hloc[i][j][1] = 0.0;
      }
    }
  }
}
    

//
//
void V_xc(fftw_ *nr, fftw_ *Vxc, double &Exc, double l_box)
{
  int i;
  int Ngr3 = Ngr*Ngr*Ngr;
  double omega = pow(l_box, 3.);
  fftw_complex *Wxc; Wxc = new fftw_complex [Ngr3];
  double asum, aasum, bsum, bbsum, exc, ndedn, rs;

  double a[4] = {0.4581652932831429, 2.217058676663745, 
                 0.7405551735357053, 0.01968227878617998};
  double b[4] = {1.0, 4.504130959426697, 1.110667363742916, 
                 0.02359291751427506};
  double dv = omega/(double) Ngr3;
  Exc = 0.0;

  for (i=0;i<Ngr3;i++) {
    rs = pow(3./4./PI/nr[i][0], 1./3.);
    asum  = a[0] +    a[1]*rs +    a[2]*rs*rs +    a[3]*rs*rs*rs;
    aasum = a[1] + 2.*a[2]*rs + 3.*a[3]*rs*rs;
    bsum  = b[0]*rs + b[1]*rs*rs + b[2]*rs*rs*rs + b[3]*pow(rs,4.0);
    bbsum = b[0] + 2.*b[1]*rs + 3.*b[2]*rs*rs + 4.*b[3]*rs*rs*rs;
    exc = -asum/bsum;
    ndedn = rs*(aasum*bsum - asum*bbsum)/bsum/bsum/3.;
    Wxc[i][0] = exc + ndedn; Wxc[i][1] = 0.0;
    Exc += exc*nr[i][0]*dv;        
  }
  fftw_plan plan;
  plan = fftw_plan_dft_3d(Ngr,Ngr,Ngr, Wxc, Vxc, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan); fftw_destroy_plan(plan); 
  
  for (i=0;i<Ngr3;i++) {
    Vxc[i][0] = Vxc[i][0]/(double) Ngr3;
    Vxc[i][1] = Vxc[i][1]/(double) Ngr3;
  }
  delete[] Wxc;
}

//
//
void H_kin(int pw[][3], fftw_ *ck, double *Hkin, double &Ekin, double l_box)
{
  int i, id, ix, iy, iz, i2;
  double G2;
  
  Ekin = 0.0;
  for (i=0;i<Npw;i++) {
    ix = pw[i][0]; iy = pw[i][1]; iz = pw[i][2];
    i2 = ix*ix + iy*iy + iz*iz;
    ix = (ix+Ngr)%Ngr; iy = (iy+Ngr)%Ngr; iz = (iz+Ngr)%Ngr;
    id = iz + Ngr*(iy + Ngr*ix);
    G2 = FOURPISQ*(double) i2 / (l_box*l_box);
    Hkin[i] = // Fill out
    Ekin += Hkin[i]*(ck[id][0]*ck[id][0] + ck[id][1]*ck[id][1]);
  }
}
