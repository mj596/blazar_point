#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "const.h"
#include "model.h"
#include "ercph.h"

/* external radiation photon field */
void
ercBEL(double r, double* ep, double* upe, int n,
       double epext, double upext, double rbel, double alpha)
{
  int i;
  
  ercBB(ep,upe,n,epext,upext);
  for (i=1;i<=n;i++)
    upe[i] *= 1.0e0/(1.0e0+(r/rbel)*(r/rbel));
    //    upe[i] *= 1.0e0/pow(r/RM,KERC);
}

void
ercIR(double r, double* ep, double* upe, int n,
      double epext, double upext, double rir)
{
  int i;
  
  ercBB(ep,upe,n,epext,upext);
  for (i=1;i<=n;i++)
    upe[i] *= 1.0e0/(1.0e0+(r/rir)*(r/rir));
    //upe[i] *= 1.0e0/((r/rir)*(r/rir));

}

void
ercPL(double* ep, double* upe, int n,
      double epext, double upext, double alpha)
{
  double ep0,epmin,epmax,normA;
  int i;
  
  ep0 = PLANCK_H*2.418e14*epext/mec2;
  normA = pow(-alpha/ep0,-alpha)*upext/ep0*1.0e0/gsl_sf_gamma(-alpha);
  epmin = (1.0e-6*ep0>PLANCK_H/mec2) ? 1.0e-6*ep0 : PLANCK_H/mec2;
  epmax = 100.0*ep0;
  for (i=1;i<=n;i++){
    ep[i]  = epmin*pow(epmax/epmin,(double)(i-1)/(double)(n-1));
    upe[i] = normA*pow(ep[i],-alpha)*exp(alpha*ep[i]/ep0);
  }
}


#define constA1 (6.52e-7)
// if frequency of maximum is estimated form vFv A1=5e-7 
#define constA2 (0.15)

void
ercBB(double* ep, double* upe, int n,
      double epext, double upext)
{
  double epmin,epmax,normA,TempX,ksi;
  int i;
  
  TempX = constA1*epext;
  normA = constA2*upext/TempX;
  epmin = (1.0e-6*TempX>PLANCK_H/mec2) ? 1.0e-6*TempX : PLANCK_H/mec2;
  epmax = 100.0*TempX;
  for (i=1;i<=n;i++){
    ep[i]  = epmin*pow(epmax/epmin,(double)(i-1)/(double)(n-1));
    ksi = ep[i]/TempX;
    upe[i] = normA*pow(ksi,3)/(exp(ksi)-1.0);
  }
}

#undef constA2
#undef constA1

/* monoenergetic version of ercph */
double
upext(double r, double up, double rbel)
{
  return up;
}
