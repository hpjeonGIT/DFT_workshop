#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"

void ck_initialize (int pw[][3], fftw_ *ck, fftw_ *nr, fftw_ *nk);
void ck_fft_cr     (fftw_ *ck, fftw_ *cr, double l_box);
void nr_fft_nk     (fftw_ *cr, fftw_ *nr, fftw_ *nk, double alpha);
void V_HT          (int pw[][3], fftw_ *nk, fftw_ **Hht, double &Eht, 
		    double l_box);
void V_loc         (int pw[][3], fftw_ *nk, fftw_ **Hloc, double &Eloc, 
		    double l_box);
void V_xc          (fftw_ *nr, fftw_ *Vxc, double &Exc, double l_box);
void H_kin         (int pw[][3], fftw_ *ck, double *Hkin, double &Ekin, 
		    double l_box);
void Hermitian_solv(int pw[][3], fftw_ **Hht, fftw_ **Hloc, fftw_ *Vxc, 
                    double *Hkin, fftw_ *ck);

int main()
{
  //
  // variable declaration
  int n;
  int pw[Npw][3] = // Fill out
  int Ngr3 = Ngr*Ngr*Ngr;
  double Ecut, Ekin, Eht, Eloc, Exc, alpha;
  double l_box = 5.0;

  double *Hkin;
  Hkin = new double [Npw];
  fftw_complex *ck, *cr, *nr, *nk, **Hht, **Hloc, *Vxc;
  ck   = new fftw_complex [Ngr3];
  cr   = new fftw_complex [Ngr3];
  nr   = new fftw_complex [Ngr3];
  nk   = new fftw_complex [Ngr3];
  Vxc  = new fftw_complex [Ngr3];
  Hht  = new fftw_complex * [Npw];
  Hloc = new fftw_complex * [Npw];
  for (n=0;n<Npw;n++) {
    Hht[n]  = new fftw_complex [Npw];
    Hloc[n] = new fftw_complex [Npw];
  }
  ck_initialize(pw, ck, nr, nk); alpha = 1.0;

  n = 0;
  do {
    ck_fft_cr(ck, cr, l_box);
    nr_fft_nk(cr, nr, nk, alpha); alpha = 0.1;
    V_HT     (pw, nk, Hht,  Eht,  l_box);
    V_loc    (pw, nk, Hloc, Eloc, l_box);
    V_xc     (    nr, Vxc,  Exc,  l_box);
    H_kin    (pw, ck, Hkin, Ekin, l_box);
    Hermitian_solv(pw, Hht, Hloc, Vxc, Hkin, ck);
    n ++;
  } while (n<100);



  std::cout << "E_kin = " << Ekin << std::endl;
  std::cout << "E_HT = " <<  Eht << std::endl;
  std::cout << "E_xc = " <<  Exc << std::endl;
  std::cout << "E_loc = " << Eloc << std::endl;


  delete[] Hkin; delete[] ck; delete[] cr;  delete[] nr; delete[] nk; 
  delete[] Vxc; 
  for (n=0;n<Npw;n++) { delete[] Hht[n]; delete[]  Hloc[n]; }
  delete[] Hht; delete[] Hloc;
  return 0;

}
