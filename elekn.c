/* solution of kinetic equation for electron number density
   synchrotron: with synchrotron self-absorption cross section
   SSC: continous energy loss approximation
   ERC: continous energy loss approximation
   Klein-Nishina included

   2005.01.10 

   2005.02.12 - any external radiation field
   2006.12.08 - Comptonization of infrared radiation added
   2009.02.16 - special functions from GNU Scientific Library

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"

#include <gsl/gsl_sf_dilog.h>

#include "const.h"
#include "model.h"
#include "ercph.h"

#define NGAMMA 200
#define NESYN  200
#define NVSYN  200
#define NVSSC  200
#define NEERC  200

#define ESP 1.0e-6

#define A43 (4.0/3.0)

double r;
double *gam,*Ng;
double *epSYN,*upeSYN;
double *epERC,*upeERC,*epERCIR,*upeERCIR;
double uB();

int
main(int argc, char* argv[])
{
  double *A,*B,*C,*S,*U;
  double epSYNMIN,epSYNMAX;
  double *vSp,vpMIN,vpMAX,*LpSvp,*vSSCp,*LpSSCvp;
  double LSv(double v, double N[]);
  double LSSCv(double v);
  double dr,err;
  double tQ(double g);
  double dgdr(double g);
  //extern void ercBEL(double r, double* epERC, double* upERC, int n, double epext, double upext, double rbel);
  //extern void ercIR(double r, double* epERCIR, double* upERCIR, int n, double epext, double upext, double rir);
  int i;
  extern void tridag(double a[],double b[],double c[],double r[],double u[],unsigned long n);
  //int input_parameters();
  extern int save(double col1[],double col2[],const char prefix[],double suffix,int n);

  double tmpN,sumN;
  FILE *output;
  double h,Linj,Lrad,eta;

  setvbuf(stdout,(char *)NULL,_IONBF,0);
  /* parameters reading */
  input_parameters(argv[1]);
  /* memory allocation */
  /* electron evolution */
  gam = dvector(0,NGAMMA+1);
  Ng = dvector(1,NGAMMA);
  A = dvector(1,NGAMMA);
  B = dvector(1,NGAMMA);
  C = dvector(1,NGAMMA);
  S = dvector(1,NGAMMA);
  U = dvector(1,NGAMMA);
  /* synchrotron radiation */
  epSYN  = dvector(1,NESYN);
  upeSYN = dvector(1,NESYN);
  vSp    = dvector(1,NVSYN);
  LpSvp  = dvector(1,NVSYN);
  /* self-synchrotron Compton */
  vSSCp   = dvector(1,NVSSC);
  LpSSCvp = dvector(1,NVSSC);
  /* external radiation field */
  epERC  = dvector(1,NEERC);
  upeERC = dvector(1,NEERC);
  epERCIR  = dvector(1,NEERC);
  upeERCIR = dvector(1,NEERC);
  /* grids initialization */
  for(i=0;i<=NGAMMA+1;i++)
    {
    //gam[i] = gammaMIN*pow(gammaMAX/gammaMIN,(double)i/(double)NGAMMA);
      gam[i] = pow(gammaMAX,(double)i/(double)NGAMMA);
      //      printf("gam[%d]=%le\n",i,gam[i]);
    }
  for(i=1;i<=NGAMMA;i++)
    A[i] = Ng[i] = 0.0;
  /* some other initializations */
  dr = (RINMAX-R0)/50.0e0;
  //  printf("dr: %le\n",dr);
  r = R0;
  epSYNMIN = PLANCK_H/mec2;  // 1Hz
  h = log(gam[2]/gam[1]);
  /* main loop */
  printf("\n\tCalculating ...\n");
  do{
    //    printf("step r=%e\n",r/R0);
    //    double KKK=uB();
    //    printf("uB: %f\n",KKK);
    if (SYN && SSC) {
      epSYNMAX = 100.0*A43*DSQR(gammaMAX)*(magB(r)/B_CR);
      for(i=1;i<=NESYN;i++){
	epSYN[i]  = epSYNMIN*pow(epSYNMAX/epSYNMIN,(i-1)/(NESYN-1.0));
	upeSYN[i] = mec2h*LSv(mec2h*epSYN[i],Ng)/(GEOMK*2.0*M_PI*(1.0-cos(THETAJ))*r*r*LIGHT_SPEED);
      }
    }
    if (ERC)
      ercBEL(r,epERC,upeERC,NEERC,GAMMA*VBEL,4.0/3.0*GAMMA*GAMMA*UBEL,RBEL,BELALPHA);
    if (ERCIR)
      ercIR(r,epERCIR,upeERCIR,NEERC,GAMMA*VIR,4.0/3.0*GAMMA*GAMMA*UIR,RIR);
    for(i=1;i<=NGAMMA;i++)
      {
	S[i] = Ng[i]+dr*tQ(gam[i]);
	//	printf("gamma[%d]=%le Ngamma[%d]=%le dr=%le tQ[%d]=%le r=%le\n",i,gam[i],i,Ng[i],dr,i,tQ(gam[i]));
      }

    //    printf("S: %e\n",S[1]);
    do {
      err = 0.0;
      /* solve equations for electrons */
      for(i=1;i<=NGAMMA;i++){
	B[i] = 1.0+dr*dgdr(0.5*(gam[i-1]+gam[i]))/(0.5*(gam[i+1]-gam[i-1]));
	C[i] = -dr*dgdr(0.5*(gam[i]+gam[i+1]))/(0.5*(gam[i+1]-gam[i-1]));
      }
      tridag(A,B,C,S,U,(unsigned long)NGAMMA);

//      int b;
//      for(b=1;b<=NGAMMA;b++)
//	{
//	  printf("B[%d]=%le C[%d]=%le S[%d]=%le U[%d]=%le\n",b,B[b],b,C[b],b,S[b],b,U[b]);
//	}

      if (SYN && SSC) {
	sumN = 0.0;
	for(i=1;i<=NESYN;i++){
	  tmpN = upeSYN[i];
	  upeSYN[i] = mec2h*LSv(mec2h*epSYN[i],U)/(GEOMK*2.0*M_PI*(1.0-cos(THETAJ))*r*r*LIGHT_SPEED);
	  sumN += upeSYN[i];
	  err += fabs(upeSYN[i] - tmpN);
	}
	err /= sumN;
	printf("iteration error %e\n",err);
      }
    } while (err>ESP);
    for(i=1;i<=NGAMMA;i++)
      Ng[i] = U[i];
    r += dr;

//    int b;
//    for(b=1;b<=NGAMMA;b++)
//      {
//	printf("Ng[%d]=%le\n",b,Ng[b]);
//      }

    save(gam,Ng,"Ng",r/R0,NGAMMA);
    /* calculate luminosities */
    double Doppler, theta;
    vpMIN = 1.0;
    if (SYN) {
      vpMAX = A43*gammaMAX*gammaMAX*(magB(r)/B_CR)*mec2h;
      vpMAX *= 100.0;
      for(i=1;i<=NVSYN;i++){
	vSp[i] = vpMIN*pow(vpMAX/vpMIN,(i-1)/(NVSYN-1.0));
	LpSvp[i] = LSv(vSp[i],Ng);
      }
      save(vSp,LpSvp,"LpS",r/R0,NVSYN);
      //      theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
      Doppler = 1.0/(GAMMA*(1.0-cos(THETAOBS)*beta(GAMMA)));
      for(i=1;i<=NVSYN;i++){
	vSp[i] *= Doppler;
	LpSvp[i] *= pow(Doppler,3);
      }
      save(vSp,LpSvp,"LS",r/R0,NVSYN);
    }
    if(SSC){
      vpMAX *= DSQR(gammaMAX);
      for(i=1;i<=NVSSC;i++){
	vSSCp[i] = vpMIN*pow(vpMAX/vpMIN,(i-1)/(NVSSC-1.0));
	LpSSCvp[i] = LSSCv(vSSCp[i]);
      }
      save(vSSCp,LpSSCvp,"LpSSC",r/R0,NVSSC);
      theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
      Doppler = 1.0/(GAMMA*(1.0-cos(THETAOBS)*beta(GAMMA)));
      printf("Doppler: %f",Doppler);
      for(i=1;i<=NVSSC;i++){
	vSSCp[i] *= Doppler;
	LpSSCvp[i] *= pow(Doppler,3);
      }
      save(vSSCp,LpSSCvp,"LSSC",r/R0,NVSSC);
    }
  }while(r<=RMAX);
  printf("\t...done.\n");
  /* free memory */
  free_dvector(upeERCIR,1,NEERC);
  free_dvector(epERCIR,1,NEERC);
  free_dvector(upeERC,1,NEERC);
  free_dvector(epERC,1,NEERC);
  free_dvector(upeSYN,1,NESYN);
  free_dvector(epSYN,1,NESYN);
  free_dvector(U,1,NGAMMA);
  free_dvector(S,1,NGAMMA);
  free_dvector(C,1,NGAMMA);
  free_dvector(B,1,NGAMMA);
  free_dvector(A,1,NGAMMA);
  free_dvector(Ng,1,NGAMMA);
  free_dvector(gam,0,NGAMMA+1);
  return 0;
}
#undef ESP
#undef NSAVE

double 
gauss(double r)
{
  return exp( -pow(r-RM,2)/(2.0*pow(RSIGMA,2)) ) ;
}

double
tri(double r)
{
  double x;
  if( r<RM )
    x = 0.5+0.5*(r-R0)/(RM-R0);
  else if( r>=RM )
    x = 1.0+0.5*(r-RM)/(RM-RINMAX);
}

/* electron injection function */ 
double
tQ(double g)
{
  double x;
  double corr = 1.0;
  
  if( QGAUSS )
    corr = gauss(r);
  
  if( QTRI )
    corr = tri(r);
  
  /* broken power law with break at gammaM */
  if(g>=gammaMIN && g<=gammaMAX && r>=R0 && r<=RINMAX)
    x = corr*ELEK*pow(g,-ELEP1)*pow(1.0e0 + pow(g/gammaM,4.0e0),(ELEP1-ELEP2)/4.0e0)/(LIGHT_SPEED*beta(GAMMA)*GAMMA);
  else
    x = 0.0e0;
  //printf("g: %le INJ: %le\n",g,x);
  //printf("%1.2e  %1.2e\n",g,x);
  return x;
}

/* electron cooling */
#define A43sigTmec (3.247899229e-8)  // 4/3 sigT/mec

double dgdr(double g)
{
  double dotg;
  double uB();
  double dotgSSC(double g);
  double dotgERC(double g);
  double dotgERCIR(double g);

  dotg = 0.0;
  /* synchrotron */
  if (SYN)
      dotg += uB();
  /* self-synchrotron-Compton */
  if (SSC)
    dotg += dotgSSC(g);
  /* external radiation Compton */
  if (ERC)
    {
      dotg += dotgERC(g);
      //      printf("dotgERC: %e\n", dotgERC(g));
    }
  if (ERCIR)
    dotg += dotgERCIR(g);
  /* multiply by common factor */
  dotg *= A43sigTmec*g*g;
  /* convert from time to radius */
  dotg *= 1.0/(beta(GAMMA)*LIGHT_SPEED*GAMMA);
  /* adiabatic cooling */
  if (ADI)
    dotg += ADIABG*g/r;
  return dotg;
}

#define A1 (3.0/8.0)
#define A2 (7.0/6.0)
#define A3 (23.0/24.0)

double
dotgSSC(double g)
     /* unnormalized FKN in total electron cooling rate due to inverse
	compton scattering of the synchrotron photon distribution
	upeSYN - integration in dlne by closed formula (4.1.14) NR */
{
  double h,b,sum;
  double fKN(double b);
  int i;
  
  /* integration step */
  h = log(epSYN[2]/epSYN[1]);
  b = 4.0*epSYN[1]*g;
  sum = A1*epSYN[1]*upeSYN[1]*fKN(b);
  b = 4.0*epSYN[2]*g;
  sum += A2*epSYN[2]*upeSYN[2]*fKN(b);
  b = 4.0*epSYN[3]*g;
  sum += A3*epSYN[3]*upeSYN[3]*fKN(b);
  for(i=4;i<=NESYN-3;i++){
    b = 4.0*epSYN[i]*g;
    sum += epSYN[i]*upeSYN[i]*fKN(b);
  }
  b = 4.0*epSYN[NESYN-2]*g;
  sum += A3*epSYN[NESYN-2]*upeSYN[NESYN-2]*fKN(b);
  b = 4.0*epSYN[NESYN-1]*g;
  sum += A2*epSYN[NESYN-1]*upeSYN[NESYN-1]*fKN(b);
  b = 4.0*epSYN[NESYN]*g;
  sum += A1*epSYN[NESYN]*upeSYN[NESYN]*fKN(b);

  return(h*sum);
}

double
dotgERC(double g)
     /* total electron cooling rate due to inverse compton scattering
	of the photon distribution upe - same as dotgSSC */
{
  double h,b,sum;
  double fKN(double b);
  int i;

  h = log(epERC[2]/epERC[1]);
  b = 4.0*epERC[1]*g;
  sum = A1*epERC[1]*upeERC[1]*fKN(b);
  b = 4.0*epERC[2]*g;
  sum += A2*epERC[2]*upeERC[2]*fKN(b);
  b = 4.0*epERC[3]*g;
  sum += A3*epERC[3]*upeERC[3]*fKN(b);
  for(i=4;i<=NEERC-3;i++){
    b = 4.0*epERC[i]*g;
    sum += epERC[i]*upeERC[i]*fKN(b);
  }
  b = 4.0*epERC[NEERC-2]*g;
  sum += A3*epERC[NEERC-2]*upeERC[NEERC-2]*fKN(b);
  b = 4.0*epERC[NEERC-1]*g;
  sum += A2*epERC[NEERC-1]*upeERC[NEERC-1]*fKN(b);
  b = 4.0*epERC[NEERC]*g;
  sum += A1*epERC[NEERC]*upeERC[NEERC]*fKN(b);
  
  return (h*sum);
}

double
dotgERCIR(double g)
/* total electron cooling rate due to inverse compton scattering
   of the photon distribution upe - same as dotgSSC */
{
  double h,b,sum;
  double fKN(double b);
  int i;

  h = log(epERCIR[2]/epERCIR[1]);
  b = 4.0*epERCIR[1]*g;
  sum = A1*epERCIR[1]*upeERCIR[1]*fKN(b);
  b = 4.0*epERCIR[2]*g;
  sum += A2*epERCIR[2]*upeERCIR[2]*fKN(b);
  b = 4.0*epERCIR[3]*g;
  sum += A3*epERCIR[3]*upeERCIR[3]*fKN(b);
  for(i=4;i<=NEERC-3;i++){
    b = 4.0*epERCIR[i]*g;
    sum += epERCIR[i]*upeERCIR[i]*fKN(b);
  }
  b = 4.0*epERCIR[NEERC-2]*g;
  sum += A3*epERCIR[NEERC-2]*upeERCIR[NEERC-2]*fKN(b);
  b = 4.0*epERCIR[NEERC-1]*g;
  sum += A2*epERCIR[NEERC-1]*upeERCIR[NEERC-1]*fKN(b);
  b = 4.0*epERCIR[NEERC]*g;
  sum += A1*epERCIR[NEERC]*upeERCIR[NEERC]*fKN(b);
  
  return (h*sum);
}

#undef A3
#undef A2
#undef A1

#define B1  (63.0/40.0)
#define B2  (441.0/200.0)
#define B3  (101.0/35.0)
#define B4  (705.0/196.0)
#define B5  (969.0/224.0)
#define B6  (3647.0/720.0)
#define B7  (1598.0/275.0)
#define B8  (3969.0/605.0)
#define B9  (8365.0/1144.0)
#define B10 (76329.0/9464.0)
#define B11 (28089.0/3185.0)
#define B12 (13403.0/1400.0)
#define B13 (112371.0/10880.0)
#define B14 (102495.0/9248.0)
#define B15 (34412.0/2907.0)

double
fKN(double b)
{
  double fg;

  if (b<0.1)
    // fg = 1.0 - B1*b + B2*b*b - B3*pow(b,3) + B4*pow(b,4) - B5*pow(b,5) + B6*pow(b,6) - B7*pow(b,7) + B8*pow(b,8) - B9*pow(b,9) + B10*pow(b,10) - B11*pow(b,11) + B12*pow(b,12) - B13*pow(b,13) + B14*pow(b,14) - B15*pow(b,15);
    fg = 1.0 - b*(B1 + b*(B2 - b*(B3 + b*(B4 - B5*b))));
  else {
    fg = ((0.5*b+6.0+6.0/b)*log(1.0+b) - ((11.0/12.0)*b*b*b + 6.0*b*b + 9.0*b + 4)/((1.0+b)*(1.0+b)) - 2.0 + 2.0e0*gsl_sf_dilog(-b));
    fg *= 9.0/(b*b*b);
  }
  return fg;
}

#undef B15
#undef B14
#undef B13
#undef B12
#undef B11
#undef B10
#undef B9
#undef B8
#undef B7
#undef B6
#undef B5
#undef B4
#undef B3
#undef B2
#undef B1

/* synchrotron luminosity */
double
LSv(double v, double N[])
{
  extern int FS(double t, double* jS, double* sigmaS);
  double a,jS,sigmaS,corr;
  double vB,h,t,LSp,tau,uB();
  double a1=3.0/8.0,a2=7.0/6.0,a3=23.0/24.0;
  int i;

  vB = ELECTRON_CHARGE*magB(r)/(2.0*M_PI*ELECTRON_MASS*LIGHT_SPEED);
  /* integration by the extended formula; eq. (4.1.14) "NR in C" */
  h = (log(gammaMAX)-log(gammaMIN))/(NGAMMA-1.0);
  LSp = tau = 0.0e0;
  t = v/(3.0*DSQR(gam[1])*vB);
  FS(t,&jS,&sigmaS);
  tau += a1*sigmaS*N[1]/pow(gam[1],4);
  LSp += a1*gam[1]*N[1]*jS;
  t = v/(3.0*DSQR(gam[2])*vB);
  FS(t,&jS,&sigmaS);
  tau += a2*sigmaS*N[2]/pow(gam[2],4);
  LSp += a2*gam[2]*N[2]*jS;
  t = v/(3.0*DSQR(gam[3])*vB);
  FS(t,&jS,&sigmaS);
  tau += a3*sigmaS*N[3]/pow(gam[3],4);
  LSp += a3*gam[3]*N[3]*jS;
  for(i=4;i<NGAMMA-2;i++){
    t = v/(3.0*DSQR(gam[i])*vB);
    FS(t,&jS,&sigmaS);
    tau += sigmaS*N[i]/pow(gam[i],4);
    LSp += gam[i]*N[i]*jS;
  }
  t = v/(3.0*DSQR(gam[NGAMMA-2])*vB);
  FS(t,&jS,&sigmaS);
  tau += a3*sigmaS*N[NGAMMA-2]/pow(gam[NGAMMA-2],4);
  LSp += a3*gam[NGAMMA-2]*N[NGAMMA-2]*jS;
  t = v/(3.0*DSQR(gam[NGAMMA-1])*vB);
  FS(t,&jS,&sigmaS);
  tau += a2*sigmaS*N[NGAMMA-1]/pow(gam[NGAMMA-1],4);
  LSp += a2*gam[NGAMMA-1]*N[NGAMMA-1]*jS;
  t = v/(3.0*DSQR(gam[NGAMMA])*vB);
  FS(t,&jS,&sigmaS);
  tau += a1*sigmaS*N[NGAMMA]/pow(gam[NGAMMA],4);
  LSp += a1*gam[NGAMMA]*N[NGAMMA]*jS;

  if (SABS && (r>R0)) {
    a = 2.0*(r-R0)/(beta(GAMMA)*GAMMA);
    tau *= h*3.0/(4.0*M_PI*a*a)*2.0*sqrt(3.0)*M_PI/15.0*ELECTRON_CHARGE/magB(r);
    corr = (tau>1.0e-5) ? (1.0-exp(-tau))/tau : (1.0-0.5*tau+1.0/6.0*tau*tau-1.0/24.0*pow(tau,3)+1.0/120.0*pow(tau,4)-1.0/720.0*pow(tau,5));
  } else
    corr = 1.0e0;
  LSp *= h*3.0*sqrt(3.0)*SIGMA_T*LIGHT_SPEED*uB()/(M_PI*vB)*corr;
  return LSp;
}

/* self synchrotron Compton luminosity */
#define constC1 (1.21e-34)  // 3/4 sigTch/mc2

double
LSSCv(double v)
{
  double dg,de,j;
  double Int1;
  int i,k;
  double fiso(double g, double e0, double e);
  double const a1=3.0/8.0,a2=7.0/6.0,a3=23.0/24.0;
  double e;

  e = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  
  dg = log(gam[2]/gam[1]);
  de = log(epSYN[2]/epSYN[1]);
  
  Int1  = a1*fiso(gam[1],epSYN[1],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[1],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[1],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[1],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[1],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[1],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[1],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  = a1*upeSYN[1]/epSYN[1]*Int1*dg;

  Int1  = a1*fiso(gam[1],epSYN[2],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[2],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[2],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[2],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[2],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[2],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[2],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  += a2*upeSYN[2]/epSYN[2]*Int1*dg;

  Int1  = a1*fiso(gam[1],epSYN[3],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[3],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[3],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[3],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[3],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[3],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[3],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  += a3*upeSYN[3]/epSYN[3]*Int1*dg;

  for(i=4;i<NESYN-2;i++) {
    Int1  = a1*fiso(gam[1],epSYN[i],e)*Ng[1]/gam[1];
    Int1 += a2*fiso(gam[2],epSYN[i],e)*Ng[2]/gam[2];
    Int1 += a3*fiso(gam[3],epSYN[i],e)*Ng[3]/gam[3];
    for(k=4;k<NGAMMA-2;k++)
      Int1 += fiso(gam[k],epSYN[i],e)*Ng[k]/gam[k];
    Int1 += a3*fiso(gam[NGAMMA-2],epSYN[i],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
    Int1 += a2*fiso(gam[NGAMMA-1],epSYN[i],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
    Int1 += a1*fiso(gam[NGAMMA],epSYN[i],e)*Ng[NGAMMA]/gam[NGAMMA];
    j += upeSYN[i]/epSYN[i]*Int1*dg;
  }

  Int1  = a1*fiso(gam[1],epSYN[NESYN-2],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[NESYN-2],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[NESYN-2],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[NESYN-2],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[NESYN-2],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[NESYN-2],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[NESYN-2],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  += a3*upeSYN[NESYN-2]/epSYN[NESYN-2]*Int1*dg;

  Int1  = a1*fiso(gam[1],epSYN[NESYN-1],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[NESYN-1],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[NESYN-1],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[NESYN-1],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[NESYN-1],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[NESYN-1],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[NESYN-1],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  += a2*upeSYN[NESYN-1]/epSYN[NESYN-1]*Int1*dg;

  Int1  = a1*fiso(gam[1],epSYN[NESYN],e)*Ng[1]/gam[1];
  Int1 += a2*fiso(gam[2],epSYN[NESYN],e)*Ng[2]/gam[2];
  Int1 += a3*fiso(gam[3],epSYN[NESYN],e)*Ng[3]/gam[3];
  for(k=4;k<NGAMMA-2;k++)
    Int1 += fiso(gam[k],epSYN[NESYN],e)*Ng[k]/gam[k];
  Int1 += a3*fiso(gam[NGAMMA-2],epSYN[NESYN],e)*Ng[NGAMMA-2]/gam[NGAMMA-2];
  Int1 += a2*fiso(gam[NGAMMA-1],epSYN[NESYN],e)*Ng[NGAMMA-1]/gam[NGAMMA-1];
  Int1 += a1*fiso(gam[NGAMMA],epSYN[NESYN],e)*Ng[NGAMMA]/gam[NGAMMA];
  j  += a1*upeSYN[NESYN]/epSYN[NESYN]*Int1*dg;
  
  j *= constC1*e*de;
  return j;
}

#undef constC1

double
fiso(double g, double e0, double e)
{
  double w,b,t,fx;
  
  if (e<g) {
    w = e/g;
    b = 4.0*e0*g;
    t = w/((1.0-w)*b);
    if((t>(1.0/(4.0*g*g))) && (t<1.0)){
      fx = 2.0*t*log(t) + t + 1.0 - 2.0*t*t + 0.5*b*t*b*t/(1.0+b*t)*(1.0-t);
      return fx;
    }
  }
  return 0.0;
}

/* magnetic field energy density */
double
uB()
{
  return(pow(magB(r),KB)/(8.0*M_PI));
}
