/* calculates lightcurve at frequency VOBS for jet source
   1999.08.31 - modified 2001.01.23
   checked 2007.01.05
   modified 2011.05.06 - single observing frequency only
   modified 2011.05.11 - sorting changed to gsl_sort and GSL interpolation used
*/

#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// #include "nrutil.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_spline.h"

#include "const.h"   
#include "model.h"

#define NTSUM 800

double lVOBS;

int
main(int argc, char* argv[])
{
  DIR *dirp;
  char *dirname;
  struct dirent *direntp;
  FILE *input,*output;
  char *filename,*t;
  float tempv,tempL;
  int i,p,q,nTIME;
  double *time;
  double *LSflare,*LSSCflare,*LERCflare,*LERCIRflare;
  double tsum[NTSUM+1],tsumMIN,tsumMAX,Lsum[NTSUM+1];
  int flare(char *name,double *time,int nTIME,double* Lvflare);

  input_parameters(argv[1]);
  lVOBS = log10(VOBS);
  /* read time */
  dirname = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(dirname,"./");
  strcat(dirname,INPUT_FILE_ID);
  strcat(dirname,"/");

  printf("dirname: %s",dirname);

  dirp = opendir(dirname);

  printf("dirname: %s",dirname);

  nTIME = 0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"LS_",(size_t)3)==0)
      nTIME++;
  time = (double *)malloc((size_t)(nTIME*sizeof(double)));
  (void)rewinddir(dirp);
  i = -1;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"LS_",(size_t)3)==0)
      time[++i] = (double)atof((direntp->d_name)+3);
  (void)closedir(dirp);
  /* sort time */
  gsl_sort(time,1,nTIME);
  
  if(SYN){
    LSflare = (double *)malloc((size_t)(nTIME*sizeof(double))) ;
    flare("LS",time,nTIME,LSflare);
  }
  if(SSC){
    LSSCflare = (double *)malloc((size_t)(nTIME*sizeof(double)));
    flare("LSSC",time,nTIME,LSSCflare);
  }
  if(ERC){
    LERCflare = (double *)malloc((size_t)(nTIME*sizeof(double)));
    flare("LERC",time,nTIME,LERCflare);
  }
  if(ERCIR){
    LERCIRflare = (double *)malloc((size_t)(nTIME*sizeof(double)));
    flare("LERCIR",time,nTIME,LERCIRflare);
  }
  /* complete flare */
  printf("\nCalculating Lflare ... ");
  fflush(stdout);

  gsl_interp_accel *acc
    = gsl_interp_accel_alloc ();
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_linear,nTIME);

  tsumMIN = time[0];
  tsumMAX = time[nTIME-1];

  for(p=0;p<=NTSUM;p++){
    tsum[p] = tsumMIN + p*(tsumMAX-tsumMIN)/(double)(NTSUM);
    Lsum[p] = 0.0e0;
  }

  if(SYN){
    gsl_spline_init(spline,time,LSflare,nTIME);
    for(p=0;p<=NTSUM;p++)
      Lsum[p] += pow(10.0e0,gsl_spline_eval(spline,tsum[p],acc));
  }

  if(SSC){
    gsl_spline_init(spline,time,LSSCflare,nTIME);
    for(p=0;p<=NTSUM;p++)
      Lsum[p] += pow(10.0e0,gsl_spline_eval(spline,tsum[p],acc));
  }

  if(ERC){
    gsl_spline_init(spline,time,LERCflare,nTIME);
    for(p=0;p<=NTSUM;p++)
      Lsum[p] += pow(10.0e0,gsl_spline_eval(spline,tsum[p],acc));
  }

  if(ERCIR){
    gsl_spline_init(spline,time,LERCIRflare,nTIME);
    for(p=0;p<=NTSUM;p++)
      Lsum[p] += pow(10.0e0,gsl_spline_eval(spline,tsum[p],acc));
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  
  char *_filename;
  _filename = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(_filename,"./");
  strcat(_filename,INPUT_FILE_ID);
  strcat(_filename,"/");
  strcat(_filename,"Lflare");
 
  output = fopen(_filename,"w");

  for(i=0;i<NTSUM;i++)
    fprintf(output,"%e %e\n",tsum[i],Lsum[i]);
  fclose(output);

  /* free memory */
  if(ERCIR)
    free(LERCIRflare);
  if(ERC)
    free(LERCflare);
  if(SSC)
    free(LSSCflare);
  if(SYN)
    free(LSflare);
  printf("done.\n");
  return 0;
}

int
flare(char *name, double *time, int nTIME, double* Lvflare)
{
  FILE *input,*output;
  char *filename,*t;
  int i,p,q,nV;
  float tempv,tempL;
  double *v,*Lv;
  
  printf("\nCalculating %sflare ... ",name);
  fflush(stdout);
  /* determine nV */
  filename = (char *)malloc((size_t)(20*sizeof(char)));
  t = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(filename,"./");
  strcat(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,name);
  strcat(filename,"_");
  sprintf(t,"%.6f",time[0]);
  strcat(filename,t);
  input = fopen(filename,"r");
  nV = 0;
  while(fscanf(input,"%f %f",&tempv,&tempL)!=EOF)
    nV++;
  fclose(input);
  /* allocate memory for v, Lv */
  v = (double *)malloc((size_t)(nV*sizeof(double)));
  Lv =(double *)malloc((size_t)(nV*sizeof(double))) ;
  /* read luminosities */
  strcpy(filename,"./");
  strcat(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,name);
  strcat(filename,"_");
  sprintf(t,"%.6f",time[0]);
  strcat(filename,t);
  input = fopen(filename,"r");
  for(p=0;p<nV;p++){
    fscanf(input,"%f %f",&tempv,&tempL);
    v[p] = (double)tempv;
  }
  fclose(input);

  /* read data and calculate flare */
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc ();
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_linear,nV);
  
  for(i=0;i<nTIME;i++){
    strcpy(filename,"./");
    strcat(filename,INPUT_FILE_ID);
    strcat(filename,"/");
    strcat(filename,name);
    strcat(filename,"_");
    sprintf(t,"%.6f",time[i]);
    strcat(filename,t);
    input = fopen(filename,"r");
    for(p=1;p<=nV;p++) {
      fscanf(input,"%f %f",&tempv,&tempL);
      Lv[p] = (double)tempL;
    }
    fclose(input);
    gsl_spline_init(spline,v,Lv,nV);
    Lvflare[i] = 0.0e0;
    if(lVOBS>v[0] && lVOBS<v[nV-1])
      Lvflare[i] = gsl_spline_eval(spline,lVOBS,acc);
  }
  free(t);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  /* save flare*/
  strcpy(filename,"./");
  strcat(filename,INPUT_FILE_ID);
  strcat(filename,"/");
  strcat(filename,name);
  strcat(filename,"flare");
  output = fopen(filename,"w");
  for(p=0;p<nTIME;p++)
    fprintf(output,"%e %e\n",time[p],Lvflare[p]);
  fclose(output);
  free(filename);

  printf("done.\n");
  return 0;
}
