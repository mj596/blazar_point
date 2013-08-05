/* calculates average spectrum for jet source
   1999.08.22 - modified 2001.01.23 */

#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#include "const.h"
#include "model.h"

#define NVSUM 300

int main(void)
{
  DIR *dirp;
  struct dirent *direntp;
  FILE *input,*output;
  int i,p,nTIME,nVS,nVSSC,nVERC,nVERCIR;
  extern void shell(unsigned long n, double a[]);
  double *time;
  double *vS,*vSSC,*vERC,*vERCIR;
  double *LSave,*LSSCave,*LERCave,*LERCIRave;
  double vsum,vsumMIN,vsumMAX,Lsum,L;
  double theta,tobsMIN;
  int average(const char *name,double *time,int nTIME,int *nV,double **v,double **Lvave);
  extern int locate(double xx[],unsigned long n,double x);
  extern double linear(double x,double x1,double y1,double x2,double y2);

  int nRAD,nG;
  double *r;
  float tempg,tempN;
  double *g,**Ng,*Ngave;
  char *filename,*rad;

  input_parameters();
  /* average electrons */
  dirp = opendir(".");
  nRAD = 0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"Ng_",(size_t)3)==0)
      nRAD++;
  if (nRAD>0){
    r = dvector(1,nRAD);
    (void)rewinddir(dirp);
    i = 0;
    while((direntp=readdir(dirp))!=NULL)
      if(strncmp(direntp->d_name,"Ng_",(size_t)3)==0)
	r[++i] = (double)atof((direntp->d_name)+3);
    shell((unsigned long)nRAD,r);
    printf("\nCalculating Ngave ... ");
    fflush(stdout);
    /* determine nG */
    filename = (char *)malloc((size_t)(20*sizeof(char)));
    rad = (char *)malloc((size_t)(20*sizeof(char)));
    strcpy(filename,"Ng_");
    sprintf(rad,"%.6f",r[1]);
    strcat(filename,rad);
    input = fopen(filename,"r");
    nG = 0;
    while(fscanf(input,"%f %f",&tempg,&tempN)!=EOF)
      nG++;
    fclose(input);
    /* allocate memory for g, Ng and Ngave */
    g = dvector(1,nG);
    Ng = dmatrix(1,nRAD,1,nG);
    Ngave = dvector(1,nG);
    /* read electron densities */
    strcpy(filename,"Ng_");
    sprintf(rad,"%.6f",r[1]);
    strcat(filename,rad);
    input = fopen(filename,"r");
    for(p=1;p<=nG;p++){
      fscanf(input,"%f %f",&tempg,&tempN);
      g[p] = (double)tempg;
      Ng[1][p] = (tempN==0.0) ? 0.0 : pow(10.0,(double)tempN);
    }
    fclose(input);
    for(i=2;i<=nRAD;i++){
      strcpy(filename,"Ng_");
      sprintf(rad,"%.6f",r[i]);
      strcat(filename,rad);
      input = fopen(filename,"r");
      for(p=1;p<=nG;p++) {
	fscanf(input,"%f %f",&tempg,&tempN);
	Ng[i][p] = (tempN==0.0) ? 0.0 : pow(10.0,(double)tempN);
      }
      fclose(input);
    }
    free(rad);  
    /* average */
    output = fopen("Ngave","w");
    for(p=1;p<=nG;p++){
      Ngave[p] = 0.0;
      for(i=1;i<=nRAD;i++)
	Ngave[p] += Ng[i][p]/nRAD;
      fprintf(output,"%e %e\n",g[p],(Ngave[p] > 0.0) ? log10(Ngave[p]) : 0.0);
    }
    fclose(output);
    free_dvector(g,1,nG);
    free_dmatrix(Ng,1,nRAD,1,nG);
    free_dvector(Ngave,1,nG);
    printf("done.\n");
  }
  (void)closedir(dirp);
  
  /* read time */
  dirp = opendir(".");
  nTIME=0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"LS_",(size_t)3)==0)
      nTIME++;
  time=dvector(1,nTIME);
  (void)rewinddir(dirp);
  i=0;
  while((direntp=readdir(dirp))!=NULL)
    if(strncmp(direntp->d_name,"LS_",(size_t)3)==0)
      time[++i]=(double)atof((direntp->d_name)+3);
  (void)closedir(dirp);
  /* sort time */
  shell((unsigned long)nTIME,time);

  /* time range */
  //  theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
  //tobsMIN = (double)R0*(1.0-beta(GAMMA)*cos(theta))/(LIGHT_SPEED*beta(GAMMA));
  /* check averaging limits */
  //if ((TAVMIN>0.0 && TAVMIN<time[1]*tobsMIN) || TAVMAX>time[nTIME]*tobsMIN)
  if ((TAVMIN>0.0 && TAVMIN<time[1]) || TAVMAX>time[nTIME])
    fprintf(stderr,"Error: averaging boundry out of range in avekn.c!\n");
  
  vsumMIN = 1.0e10; /* default values */
  vsumMAX = 1.0e20;
  if(SYN){
    average("LS",time,nTIME,&nVS,&vS,&LSave);
    if (vS[1]<vsumMIN)
      vsumMIN = vS[1];
    if (vS[nVS]>vsumMAX)
      vsumMAX = vS[nVS];
  }
  if(SSC){
    average("LSSC",time,nTIME,&nVSSC,&vSSC,&LSSCave);
    if (vSSC[1]<vsumMIN)
      vsumMIN = vSSC[1];
    if (vSSC[nVSSC]>vsumMAX)
      vsumMAX = vSSC[nVSSC];
  }
  if(ERC){
    average("LERC",time,nTIME,&nVERC,&vERC,&LERCave);
    if (vERC[1]<vsumMIN)
      vsumMIN = vERC[1];
    if (vERC[nVERC]>vsumMAX)
      vsumMAX = vERC[nVERC];
  }
  if(ERCIR){
    average("LERCIR",time,nTIME,&nVERCIR,&vERCIR,&LERCIRave);
    if (vERCIR[1]<vsumMIN)
      vsumMIN = vERCIR[1];
    if (vERCIR[nVERCIR]>vsumMAX)
      vsumMAX = vERCIR[nVERCIR];
  }
  /* complete averaged spectra */
  output = fopen("Lave","w");
  printf("\nCalculating Lave ... ");
  fflush(stdout);
  for(i=1;i<=NVSUM;i++){
    vsum = vsumMIN*pow(vsumMAX/vsumMIN,(i-1.0)/(NVSUM-1.0));
    Lsum = 0.0;
    if(SYN && vsum>vS[1] && vsum<vS[nVS]){
      p = locate(vS,nVS,vsum);
      L = linear(vsum,vS[p],LSave[p],vS[p+1],LSave[p+1]);
      Lsum += L;
    }
    if(SSC && vsum>vSSC[1] && vsum<vSSC[nVSSC]){
      p = locate(vSSC,nVSSC,vsum);
      L = linear(vsum,vSSC[p],LSSCave[p],vSSC[p+1],LSSCave[p+1]);
      Lsum += L;
    }
    if(ERC && vsum>vERC[1] && vsum<vERC[nVERC]){
      p = locate(vERC,nVERC,vsum);
      L = linear(vsum,vERC[p],LERCave[p],vERC[p+1],LERCave[p+1]);
      Lsum += L;
    }
    if(ERCIR && vsum>vERCIR[1] && vsum<vERCIR[nVERCIR]){
      p = locate(vERCIR,nVERCIR,vsum);
      L = linear(vsum,vERCIR[p],LERCIRave[p],vERCIR[p+1],LERCIRave[p+1]);
      Lsum += L;
    }
    fprintf(output,"%e %e\n",log10(vsum),(Lsum==0.0) ? 0.0 : log10(Lsum));
  }
  fclose(output);
  printf("done.\n");
  /* free memory */
  if(ERCIR){
    free_dvector(LERCIRave,1,nVERCIR);
    free_dvector(vERCIR,1,nVERCIR);
  }
  if(ERC){
    free_dvector(LERCave,1,nVERC);
    free_dvector(vERC,1,nVERC);
  }
  if(SSC){
    free_dvector(LSSCave,1,nVSSC);
    free_dvector(vSSC,1,nVSSC);
  }
  if(SYN){
    free_dvector(LSave,1,nVS);
    free_dvector(vS,1,nVS);
  }
  free_dvector(time,1,nTIME);
  return 0;
}

int
average(const char *name, double *time, int nTIME, int *nV, double **v, double **Lvave)
{
  FILE *input,*output;
  char *filename,*t;
  int i,p,t1,t2;
  float tempv,tempL;
  double **Lv,dt;
  extern int locate(double xx[],unsigned long n,double x);

  printf("\nCalculating %save ... ",name);
  fflush(stdout);
  /* determine nV */
  filename=(char *)malloc((size_t)(20*sizeof(char)));
  t=(char *)malloc((size_t)(20*sizeof(char)));
  strcpy(filename,name);
  strcat(filename,"_");
  sprintf(t,"%.6f",time[1]);
  strcat(filename,t);
  input = fopen(filename,"r");
  *nV=0;
  while(fscanf(input,"%f %f",&tempv,&tempL)!=EOF)
    (*nV)++;
  fclose(input);
  /* allocate memory for v, Lv and Lvave */
  (*v)=dvector(1,(*nV));
  Lv=dmatrix(1,nTIME,1,(*nV));
  (*Lvave)=dvector(1,(*nV));
  /* read luminosities */
  strcpy(filename,name);
  strcat(filename,"_");
  sprintf(t,"%.6f",time[1]);
  strcat(filename,t);
  input=fopen(filename,"r");
  for(p=1;p<=(*nV);p++){
    fscanf(input,"%f %f",&tempv,&tempL);
    (*v)[p]=pow(10.0,(double)tempv);
    Lv[1][p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
  }
  fclose(input);
  for(i=2;i<=nTIME;i++){
    strcpy(filename,name);
    strcat(filename,"_");
    sprintf(t,"%.6f",time[i]);
    strcat(filename,t);
    input=fopen(filename,"r");
    for(p=1;p<=(*nV);p++) {
      fscanf(input,"%f %f",&tempv,&tempL);
      Lv[i][p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
    }
    fclose(input);
  }
  free(t);
  /* do average and save */
  strcpy(filename,name);
  strcat(filename,"ave");
  /* averaging boundry */
  if (TAVMIN<0.0)
    t1 = 1;
  else {
    t1 = locate(time,nTIME,TAVMIN);
    if (t1==0 || t1==nTIME)
      t1 = 1;
  }
  if (TAVMAX<0.0)
    t2 = nTIME;
  else {
    t2 = locate(time,nTIME,TAVMAX);
    if (t2==0 || t2==nTIME)
      t2 = nTIME;
  }
  dt = t2-t1+1.0;
  output=fopen(filename,"w");
  for(p=1;p<=(*nV);p++){
    (*Lvave)[p]=0.0;
    for(i=t1;i<=t2;i++)
      (*Lvave)[p] += Lv[i][p]/dt;
    fprintf(output,"%e %e\n",log10((*v)[p]),((*Lvave)[p]==0.0) ? 0.0 : log10((*Lvave)[p]));
  }
  fclose(output);
  /* free memory */
  free_dmatrix(Lv,1,nTIME,1,(*nV));
  free(filename);
  printf("done.\n");
  return 0;
}
