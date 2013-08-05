/* calculates external compton luminosity for point source
   2004.09.29 */

#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#include "const.h"
#include "model.h"

#define NVERC 200

#define NEERC 200
double *epERC,*upeERC;

int
main(int argc, char* argv[])
{
  DIR *dirp;
  char *dirname;
  struct dirent *direntp;  
  FILE *input;
  char *filename,*rad;
  float tempg,tempN;
  int i,k,p;
  int nRAD,nGAMMA;
  double *r,*gamma,**Ng;
  double *vp,vpMIN,vpMAX,*Lvp;
  extern void shell(unsigned long int n, double a[]);
  int input_parameters();
  extern int save(double col1[],double col2[],const char prefix[],double suffix,int n);
  double LERCv(double const v,double const g[],double const N[],const int n,double const theta);
  extern void ercph(double r, double* epERC, double* upERC, int n, double epext, double upext, double rbel);
  
  /* read parameters */
  input_parameters(argv[1]);
  /* read radii */
  dirname = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(dirname,"./");
  strcat(dirname,INPUT_FILE_ID);
  strcat(dirname,"/");
  dirp = opendir(dirname);
  nRAD = 0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"Ng_",(size_t)3) == 0)
      nRAD++;
  r = dvector(1,nRAD);
  (void)rewinddir(dirp);
  i = 0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"Ng_",(size_t)3) == 0)
      r[++i] = (double)R0*atof((direntp->d_name)+3);
  (void)closedir(dirp);
  /* sort radii */
  shell((unsigned long)nRAD,r);

  for(i=1;i<nRAD;i++)
    printf("rad[%d]: %le\n",i,r[i]);

  if(ERC) {
    /* read gamma & N_gamma */
    printf("\nReading Ng ... ");
    fflush(stdout);
    /* determine nGAMMA */
    filename = (char *)malloc((size_t)(20*sizeof(char)));
    rad = (char *)malloc((size_t)(20*sizeof(char)));
    strcpy(filename,"./");
    strcat(filename,INPUT_FILE_ID);
    strcat(filename,"/");
    strcat(filename,"Ng");
    strcat(filename,"_");
    sprintf(rad,"%.6f",r[1]/R0);
    strcat(filename,rad);
    input = fopen(filename,"r");
    nGAMMA = 0;
    while(fscanf(input,"%f %f",&tempg,&tempN)!=EOF)
      nGAMMA++;
    fclose(input);
    printf("nGamma: %d\n",nGAMMA);
    /* allocate memory for gamma and Ng */
    gamma = dvector(1,nGAMMA);
    Ng = dmatrix(1,nRAD,1,nGAMMA);
    /* read gamma and Ng */
    strcpy(filename,"./");
    strcat(filename,INPUT_FILE_ID);
    strcat(filename,"/");
    strcat(filename,"Ng");
    strcat(filename,"_");
    sprintf(rad,"%.6f",r[1]/R0);
    strcat(filename,rad);
    input = fopen(filename,"r");
    printf("filename: %s\n",filename);
    for(p=1;p<=nGAMMA;p++) {
      fscanf(input,"%f %f",&tempg,&tempN);
      gamma[p] = pow(10.0,(double)tempg);
      Ng[1][p] = (tempN==0.0) ? 0.0 : pow(10.0,(double)tempN);
    }
    fclose(input);
    for(i=2;i<=nRAD;i++){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"Ng");
      strcat(filename,"_");
      sprintf(rad,"%.6f",r[i]/R0);
      strcat(filename,rad);
      input = fopen(filename,"r");
      printf("filename: %s\n",filename);
      for(p=1;p<=nGAMMA;p++) {
	fscanf(input,"%f %f",&tempg,&tempN);
	Ng[i][p] = (tempN==0.0) ? 0.0 : pow(10.0,(double)tempN);
      }
      fclose(input);
    }
    free(rad);
    free(filename);
    printf("done.\n");
    double Doppler, theta;

    /* calculate ERC luminosity */
    printf("Calculating LERC ... ");
    fflush(stdout);
    /* external field */
    epERC = dvector(1,NEERC);
    upeERC = dvector(1,NEERC);
    /* allocate memory */
    vp = dvector(1,NVERC);
    Lvp = dvector(1,NVERC);
    for(k=1;k<=nRAD;k++){
      ercBEL(r[k],epERC,upeERC,NEERC,GAMMA*VBEL,4.0/3.0*GAMMA*GAMMA*UBEL,RBEL,BELALPHA);
      printf("epERC[1]: %le\n",epERC[1]);
      printf("epERC[NEERC]: %le\n",epERC[NEERC]);
      vpMIN = 0.9*mec2h*epERC[1];
      vpMAX = 1.1*DSQR(gamma[nGAMMA])*mec2h*epERC[NEERC];
      for(i=1;i<=NVERC;i++){
	vp[i] = vpMIN*pow(vpMAX/vpMIN,(i-1.0)/(NVERC-1.0));
	Lvp[i] = LERCv(vp[i],gamma,Ng[k],nGAMMA,THETAOBS);
//	printf("vp: %le\n",vp[i]);
//	printf("gamma: %le\n",gamma);
//	printf("Ng: %le\n",Ng[k]);
//	printf("THETAOBS: %le\n",THETAOBS);
//	printf("nGAMMA: %le\n",GAMMA);
//	printf("Lvp: %le\n",Lvp[i]);
      }

      //   theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
      Doppler = 1.0/(GAMMA*(1.0-cos(THETAOBS)*beta(GAMMA)));
      printf("Doppler: %le\n",Doppler);
      for(i=1;i<=NVERC;i++){
	vp[i] *= Doppler;
	Lvp[i] *= pow(Doppler,3);
      }
      save(vp,Lvp,"LERC",r[k]/R0,NVERC);
    }
    /* free memory */
    free_dvector(upeERC,1,NEERC);
    free_dvector(epERC,1,NEERC);
    free_dvector(vp,1,NVERC);
    free_dvector(Lvp,1,NVERC);
    printf("done\n");
    free_dmatrix(Ng,1,nRAD,1,nGAMMA);
    free_dvector(gamma,1,nGAMMA);
  }
  free_dvector(r,1,nRAD);
  return 0;
}

/* external Compton luminosity */
#define constC1 (1.21e-34)  // 3/4 sigTch/mc2

double
LERCv(double const v, double const g[], double const N[], int const n, double const theta)
{
  double dg,de,j;
  double Int1;
  int i,k;
  double f(double g, double e0, double e, double miu);
  double const a1=3.0/8.0,a2=7.0/6.0,a3=23.0/24.0;
  double e,miu;

  e = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  miu = -(cos(theta)-beta(GAMMA))/(1.0-beta(GAMMA)*cos(theta));

  //int z;
  //  for(z=1;z<n;z++)
  //  printf("gamma[%d]=%le\n",z,g[z]);

  dg = log(g[2]/g[1]);
  de = log(epERC[2]/epERC[1]);
  Int1  = a1*f(g[1],epERC[1],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[1],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[1],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[1],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[1],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[1],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[1],e,miu)*N[n]/g[n];
  j  = a1*upeERC[1]/epERC[1]*Int1*dg;

  Int1  = a1*f(g[1],epERC[2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[2],e,miu)*N[n]/g[n];
  j  += a2*upeERC[2]/epERC[2]*Int1*dg;

  Int1  = a1*f(g[1],epERC[3],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[3],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[3],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[3],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[3],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[3],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[3],e,miu)*N[n]/g[n];
  j  += a3*upeERC[3]/epERC[3]*Int1*dg;

  for(i=4;i<NEERC-2;i++) {
    Int1  = a1*f(g[1],epERC[i],e,miu)*N[1]/g[1];
    Int1 += a2*f(g[2],epERC[i],e,miu)*N[2]/g[2];
    Int1 += a3*f(g[3],epERC[i],e,miu)*N[3]/g[3];
    for(k=4;k<n-2;k++)
      Int1 += f(g[k],epERC[i],e,miu)*N[k]/g[k];
    Int1 += a3*f(g[n-2],epERC[i],e,miu)*N[n-2]/g[n-2];
    Int1 += a2*f(g[n-1],epERC[i],e,miu)*N[n-1]/g[n-1];
    Int1 += a1*f(g[n],epERC[i],e,miu)*N[n]/g[n];
    j += upeERC[i]/epERC[i]*Int1*dg;
  }

  Int1  = a1*f(g[1],epERC[NEERC-2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC-2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC-2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC-2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC-2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC-2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC-2],e,miu)*N[n]/g[n];
  j  += a3*upeERC[NEERC-2]/epERC[NEERC-2]*Int1*dg;

  Int1  = a1*f(g[1],epERC[NEERC-1],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC-1],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC-1],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC-1],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC-1],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC-1],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC-1],e,miu)*N[n]/g[n];
  j  += a2*upeERC[NEERC-1]/epERC[NEERC-1]*Int1*dg;

  Int1  = a1*f(g[1],epERC[NEERC],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC],e,miu)*N[n]/g[n];
  j  += a1*upeERC[NEERC]/epERC[NEERC]*Int1*dg;
  
  j *= constC1*e*de;
  return j;
}

#undef contC1

double
f(double g, double e0, double e, double miu)
{
  double w,wp,b,t,fx;
  
  w = e/g;
  b = 2.0*(1.0-miu)*e0*g;
  if (e>e0 && e<b*g/(1+b)) {
    //    printf("gamma: %le e0: %le e: %le miu:  %le\n",g,e0,e,miu);
    wp = 1.0-w;
    fx = (1.0+0.5*w*w/wp-2.0*w/(b*wp)+2.0*w*w/(b*b*wp*wp));
    if (fx>0.0)
      return fx;
    else
      return 0.0;
  } else
    return 0.0;
}
