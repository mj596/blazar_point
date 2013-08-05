#include <gsl/gsl_sf_bessel.h>

int
FS(double t, double* jS, double* sigmaS)
{
  double K13,K43;

  if(t<50.0){
    K13 = gsl_sf_bessel_Knu(1.0e0/3.0e0,t);
    K43 = gsl_sf_bessel_Knu(4.0e0/3.0e0,t);
    (*sigmaS) = t*(K43*K43-K13*K13);
    (*jS) = t*t*(K43*K13-3.0e0/5.0e0*(*sigmaS));
  } else
    (*sigmaS) = (*jS) = 0.0;
  return 0;
}
