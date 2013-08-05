#include "nrutil.h"

void tridag(double a[], double b[], double c[], double r[], 
	    double u[], unsigned long n)
  /* Solves for a vector u[1..n] of length n the tridiagonal linear set 
     of equations. a[1..n], b[1..n], c[1..n], and r[1..n] are input 
     vectors */
{
  unsigned long i;
  double bet,*gam;
  
  gam=dvector(1,n);
  if(b[1]==0.0)
    nrerror("tridag: rewrite equations");
  u[1]=r[1]/(bet=b[1]);
  for(i=2;i<=n;i++){
    gam[i]=c[i-1]/bet;
    bet=b[i]-a[i]*gam[i];
    if(bet==0.0)
      nrerror("tridag: crash");
    u[i]=(r[i]-a[i]*u[i-1])/bet;
  }
  for(i=(n-1);i>=1;i--)
    u[i]-=gam[i+1]*u[i+1];
  free_dvector(gam,1,n);
}
