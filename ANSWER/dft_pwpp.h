/*
Basic Unit
1 Hartree              = 27.211 383 86(68) eV   () means uncertainty
1 bohr radius a0       = 0.529 177 208 59(36) x 10-10 m    
1 electron mass (me)   = 9.109 382 15(45) x 10-31 kg    
1 elementary charge    = 1.602 176 487(40) x 10-19 C    
electric constant eps0 = 8.854 187 817... x 10-12 F m-1    
planck constant h      = 4.135 667 33(10) x 10-15 eV s    
h/2pi  hbar            = 1.054 571 628(53) x 10-34 J s 
*/
#define PI 3.1415926535897931
#define TWOPI 6.2831853071795862     // 2.*PI
#define FOURPISQ 39.478417604357432  // (2PI)^2
#define Npw 7
#define Ngr 5
#include <vector>

typedef struct{
  double real, imag;
} dcomplex;

typedef fftw_complex fftw_;
