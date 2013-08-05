/* calculates luminosities for jet source
   full Klein-Nishina version
   2005.02.01 */

#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#include "const.h"
#include "model.h"
#include "ercph.h"

#define NMI 50

#define NVS   200
#define NVSSC 200
#define NVERC 200

#define NEERC 200
double *epERC,*upeERC,*epERCIR,*upeERCIR;

#define NTIME 10

int
main(int argc, char* argv[])
{
  DIR *dirp;
  char *dirname;
  struct dirent *direntp;  
  FILE *input;
  char *filename,*rad;
  float tempg,tempN,tempv,tempL;
  int i,k,l,p,q,q1,q2;
  extern int locate(double xx[],unsigned long n,double x);
  extern double linear(double x,double x1,double y1,double x2,double y2);
  int nRAD,nVSP,nVSSCP,nGAMMA;
  double *r,robs,*tobs,tobsMIN,tobsMAX,t,u,findr(double t, double theta);
  double *gamma,**Ng,*N;
  double vpobs,*vobs,vobsMIN,vobsMAX,**vSp,**vSSCp;
  double Lvp,*vLvobs,vLvp,*Lvobs,**LpSvp,**LpSSCvp;
  double theta,dtheta,dfi;
  double gg,Ngg,Doppler;
  extern void shell(unsigned long int n, double a[]);
  int input_parameters();
  extern int save(double col1[],double col2[],const char prefix[],double r,int n);
  double *vpERC,*vpERCIR,**LpERCvp,**LpERCvpIR;
  double LERCv(double const v,double const g[],double const N[],const int n,double const theta);
  double LERCIRv(double const v,double const g[],double const N[],const int n,double const theta);

  setvbuf(stdout,(char *)NULL,_IONBF,0);
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

  /* time range */
  theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
  tobsMIN = (double)R0*(1.0-beta(GAMMA)*cos(theta))/(LIGHT_SPEED*beta(GAMMA));
  tobsMAX = r[nRAD]*(1.0-beta(GAMMA)*cos(theta))/(LIGHT_SPEED*beta(GAMMA));
  tobs = dvector(0,NTIME);
  for(i=0;i<=NTIME;i++)
    tobs[i] = tobsMIN + (double)i*(tobsMAX-tobsMIN)/(double)NTIME;
  printf("\n\ttobs0 = %e\n",tobs[0]);

  /* calculations */
  if(SYN){
    /* calculate SYN luminosity */
    printf("\nReading LpS ... ");
    fflush(stdout);
    /* determine nVSP */
    filename = (char *)malloc((size_t)(20*sizeof(char)));
    rad = (char *)malloc((size_t)(20*sizeof(char)));
    strcpy(filename,"./");
    strcat(filename,INPUT_FILE_ID);
    strcat(filename,"/");
    strcat(filename,"LpS");
    strcat(filename,"_");
    sprintf(rad,"%.6f",r[1]/R0);
    strcat(filename,rad);
    input = fopen(filename,"r");
    nVSP = 0;
    while(fscanf(input,"%f %f",&tempv,&tempL)!=EOF)
      nVSP++;  
    fclose(input);
    /* allocate memory for vSp and LpSvp */
    vSp = dmatrix(1,nRAD,1,nVSP);
    LpSvp = dmatrix(1,nRAD,1,nVSP);
    /* read internal synchrotron luminosities */
    for(i=1;i<=nRAD;i++){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LpS");
      strcat(filename,"_");
      sprintf(rad,"%.6f",r[i]/R0);
      strcat(filename,rad);
      input = fopen(filename,"r");
      for(p=1;p<=nVSP;p++) {
        fscanf(input,"%f %f",&tempv,&tempL);
        vSp[i][p] = pow(10.0,(double)tempv);
        LpSvp[i][p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
    }
    free(rad);
    free(filename);
    printf("done.\n");
    
    printf("Calculating LS ... ");
    fflush(stdout);
    vobsMIN = 1.0e30;
    vobsMAX = 0.0;
    vobs = dvector(1,NVS);
    Lvobs = dvector(1,NVS);
    theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
    for(i=1;i<=nRAD;i++)
      if(vobsMAX<Doppler*vSp[i][nVSP])
        vobsMAX = Doppler*vSp[i][nVSP];
    theta = THETAOBS+THETAJ;
    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
    for(i=1;i<=nRAD;i++)
      if(vobsMIN>Doppler*vSp[i][1])
        vobsMIN = Doppler*vSp[i][1];
    /* set frequency grid */
    for(i=1;i<=NVS;i++)
      vobs[i] = vobsMIN*pow(vobsMAX/vobsMIN,(i-1.0)/(NVS-1.0));
    dtheta = 2.0e0*THETAJ/NMI;

    for (k=0;k<=NTIME;k++)
      {
	for(i=1;i<=NVS;i++)
	  {
	    Lvobs[i] = 0.0e0;
	    theta = THETAOBS + THETAJ;      
	    do
	      {
		Lvp = 0.0e0;
		robs = findr(tobs[k],theta);
          
		if(robs>r[1] && robs<r[nRAD])
		  {
		    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
		    vpobs = vobs[i]/Doppler;
		    p = locate(r,nRAD,robs);

		    if(vpobs>DMAX(vSp[p][1],vSp[p+1][1]) && vpobs<DMIN(vSp[p][nVSP],vSp[p+1][nVSP]))
		      {
			/* 2D interpolation - modified bilinear */
			q1 = locate(vSp[p],nVSP,vpobs);
			q2 = locate(vSp[p+1],nVSP,vpobs);
			t = sqrt((vpobs-vSp[p][q1])*(vpobs-vSp[p+1][q2])/((vSp[p][q1+1]-vSp[p][q1])*(vSp[p+1][q2+1]-vSp[p+1][q2])));
			u = (robs-r[p])/(r[p+1]-r[p]);
			Lvp = (1.0-t)*(1.0-u)*LpSvp[p][q1]+t*(1.0-u)*LpSvp[p][q1+1]+t*u*LpSvp[p+1][q2+1]+(1.0-t)*u*LpSvp[p+1][q2];
			if(theta<(THETAJ-THETAOBS) || THETAOBS==0.0)
			  dfi = M_PI;
			else
			  {
			    if (fabs(((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS))))>1.0e0)
			      dfi = 2.0e0*M_PI;
			    else
			      dfi = 2.0e0*acos((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS)));
			  }
			Lvobs[i] += dfi*Lvp*pow(Doppler,3)*fabs(sin(theta));
		      }
		  }
		theta -= dtheta;
	      } while(theta>(THETAOBS-THETAJ));
	    Lvobs[i] *= dtheta/(2.0*M_PI*(1.0-cos(THETAJ)));
	  }

      save(vobs,Lvobs,"LS",tobs[k]/tobs[0],NVS);
    }
    /* free memory */
    free_dvector(vobs,1,NVS);
    free_dvector(Lvobs,1,NVS);
    free_dmatrix(LpSvp,1,nRAD,1,nVSP);
    free_dmatrix(vSp,1,nRAD,1,nVSP);
    printf("done\n");
  }
  
  if(SSC){
    /* calculate SSC luminosity */
    printf("\nReading LpSSC ... ");
    fflush(stdout);
    /* determine nVSSCP */  
    filename = (char *)malloc((size_t)(20*sizeof(char)));
    rad = (char *)malloc((size_t)(20*sizeof(char)));
    strcpy(filename,"./");
    strcat(filename,INPUT_FILE_ID);
    strcat(filename,"/");
    strcat(filename,"LpSSC");
    strcat(filename,"_");
    sprintf(rad,"%.6f",r[1]/R0);
    strcat(filename,rad);
    input = fopen(filename,"r");
    nVSSCP = 0;
    while(fscanf(input,"%f %f",&tempv,&tempL)!=EOF)
      nVSSCP++;  
    fclose(input);
    /* allocate memory for vSSCp and LpSSCvp */
    vSSCp = dmatrix(1,nRAD,1,nVSSCP);
    LpSSCvp = dmatrix(1,nRAD,1,nVSSCP);
    /* read internal SSC luminosities */
    for(i=1;i<=nRAD;i++){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LpSSC");
      strcat(filename,"_");
      sprintf(rad,"%.6f",r[i]/R0);
      strcat(filename,rad);
      input = fopen(filename,"r");
      for(p=1;p<=nVSSCP;p++) {
	fscanf(input,"%f %f",&tempv,&tempL);
	vSSCp[i][p] = pow(10.0,(double)tempv);
	LpSSCvp[i][p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
    }
    free(rad);
    free(filename);
    printf("done.\n");

    printf("Calculating LSSC ... ");
    fflush(stdout);
    vobs = dvector(1,NVSSC);
    Lvobs = dvector(1,NVSSC);
    vobsMIN = 1.0e33;
    vobsMAX = 0.0;
    vobs = dvector(1,NVS);
    Lvobs = dvector(1,NVS);
    theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
    for(p=1;p<=nRAD;p++)
      if(vobsMAX<Doppler*vSSCp[p][nVSSCP])
	vobsMAX = Doppler*vSSCp[p][nVSSCP];
    theta = THETAOBS+THETAJ;
    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
    for(p=1;p<=nRAD;p++)
      if(vobsMIN>Doppler*vSSCp[p][1])
	vobsMIN = Doppler*vSSCp[p][1];
    for(i=1;i<=NVSSC;i++)
      vobs[i] = vobsMIN*pow(vobsMAX/vobsMIN,(i-1.0)/(NVSSC-1.0));
    dtheta = 2.0*THETAJ/NMI;
    for (k=0;k<=NTIME;k++){
      for(i=1;i<=NVSSC;i++){
	Lvobs[i] = 0.0;
	theta = THETAOBS+THETAJ;      
	do{
	  robs = findr(tobs[k],theta);
	  if(robs>r[1] && robs<r[nRAD]){
	    Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
	    vpobs = vobs[i]/Doppler;
	    p = locate(r,nRAD,robs);
	    if(vpobs>DMAX(vSSCp[p][1],vSSCp[p+1][1]) && vpobs<DMIN(vSSCp[p][nVSSCP],vSSCp[p+1][nVSSCP])){
	      /* 2D interpolation - modified bilinear */
	      q1 = locate(vSSCp[p],nVSSCP,vpobs);
	      q2 = locate(vSSCp[p+1],nVSSCP,vpobs);
	      t = sqrt((vpobs-vSSCp[p][q1])*(vpobs-vSSCp[p+1][q2])/((vSSCp[p][q1+1]-vSSCp[p][q1])*(vSSCp[p+1][q2+1]-vSSCp[p+1][q2])));
	      u = (robs-r[p])/(r[p+1]-r[p]);
	      Lvp = (1.0-t)*(1.0-u)*LpSSCvp[p][q1]+t*(1.0-u)*LpSSCvp[p][q1+1]+t*u*LpSSCvp[p+1][q2+1]+(1.0-t)*u*LpSSCvp[p+1][q2];
	      if(theta<(THETAJ-THETAOBS) || THETAOBS==0.0)
		dfi = M_PI;
	      else {
		if (fabs(((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS))))>1.0e0)
		  dfi = 2.0e0*M_PI;
		else
		  dfi = 2.0e0*acos((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS)));
	      }
	      Lvobs[i] += dfi*Lvp*pow(Doppler,3)*fabs(sin(theta));
	    }
	  }
	  theta -= dtheta;
	}while(theta>=(THETAOBS-THETAJ));
	Lvobs[i] *= dtheta/(2.0*M_PI*(1.0-cos(THETAJ)));
      }
      save(vobs,Lvobs,"LSSC",tobs[k]/tobs[0],NVSSC);
    }
    /* free memory */
    free_dvector(vobs,1,NVSSC);
    free_dvector(Lvobs,1,NVSSC);
    free_dmatrix(LpSSCvp,1,nRAD,1,nVSSCP);
    free_dmatrix(vSSCp,1,nRAD,1,nVSSCP);
    printf("done\n");
  }

  if(ERC || ERCIR){
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
    for(p=1;p<=nGAMMA;p++) {
      fscanf(input,"%f %f",&tempg,&tempN);
      gamma[p]=pow(10.0,(double)tempg);
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
      input=fopen(filename,"r");
      for(p=1;p<=nGAMMA;p++) {
	fscanf(input,"%f %f",&tempg,&tempN);
	Ng[i][p] = (tempN==0.0) ? 0.0 : pow(10.0,(double)tempN);
      }
      fclose(input);
    }
    free(rad);
    free(filename);
    printf("done.\n");

    /* external field - broad emission lines */
    if (ERC) {
      epERC = dvector(1,NEERC);
      upeERC = dvector(1,NEERC);
      ercBEL(r[1],epERC,upeERC,NEERC,GAMMA*VBEL,4.0/3.0*GAMMA*GAMMA*UBEL,RBEL,BELALPHA);
      /* calculate ERC luminosity */
      printf("Calculating LERC ... ");
      fflush(stdout);
      /* allocate memory */
      vobs = dvector(1,NVERC);
      Lvobs = dvector(1,NVERC);
      theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
      vobsMAX = 10.0*DSQR(Doppler*gamma[nGAMMA])*epERC[NEERC]*mec2h;
      if (vobsMAX>100.0*gamma[nGAMMA]*mec2h)
	vobsMAX = 100.0*gamma[nGAMMA]*mec2h;
      theta = THETAOBS+THETAJ;
      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
      vobsMIN = DSQR(Doppler*gamma[1])*epERC[1]*mec2h;
      for(i=1;i<=NVERC;i++)
	vobs[i]=vobsMIN*pow(vobsMAX/vobsMIN,(i-1.0)/(NVERC-1.0));
      N = dvector(1,nGAMMA);
      dtheta = 2.0*THETAJ/NMI;
      for(k=0;k<=NTIME;k++){
	for(i=1;i<=NVERC;i++){
	  Lvobs[i] = 0.0;
	  theta = THETAOBS+THETAJ;      
	  do{
	    robs = findr(tobs[k],theta);
	    if(robs>r[1] && robs<r[nRAD]){
	      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
	      vpobs = vobs[i]/Doppler;
	      p = locate(r,nRAD,robs);
	      // interpolate Ng into N
	      for(l=1;l<=nGAMMA;l++)
		if (Ng[p][l]>0.0 && Ng[p+1][l]>0.0)
		  N[l] = exp(linear(log(robs),log(r[p]),log(Ng[p][l]),log(r[p+1]),log(Ng[p+1][l])));
		else
		  N[l] = 0.0;
	      ercBEL(robs,epERC,upeERC,NEERC,GAMMA*VBEL,4.0/3.0*GAMMA*GAMMA*UBEL,RBEL,BELALPHA);
	      Lvp = LERCv(vpobs,gamma,N,nGAMMA,theta);
	      if(theta<(THETAJ-THETAOBS) || THETAOBS==0.0)
		dfi = M_PI;
	      else {
		if (fabs(((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS))))>1.0e0)
		  dfi = 2.0e0*M_PI;
		else
		  dfi = 2.0e0*acos((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS)));
	      }
	      Lvobs[i] += dfi*Lvp*pow(Doppler,3)*fabs(sin(theta));
	    }
	    theta -= dtheta;
	  }while(theta>=(THETAOBS-THETAJ));
	  Lvobs[i] *= dtheta/(2.0*M_PI*(1.0-cos(THETAJ)));
	}
	save(vobs,Lvobs,"LERC",tobs[k]/tobs[0],NVERC);
      }
      /* free memory */
      free_dvector(N,1,nGAMMA);
      free_dvector(Lvobs,1,NVERC);
      free_dvector(vobs,1,NVERC);
      free_dvector(upeERC,1,NEERC);
      free_dvector(epERC,1,NEERC);
      printf("done\n");
    }

    /* external field - infrared */
    if (ERCIR) {
      epERCIR = dvector(1,NEERC);
      upeERCIR = dvector(1,NEERC);
      ercIR(r[1],epERCIR,upeERCIR,NEERC,GAMMA*VIR,4.0/3.0*GAMMA*GAMMA*UIR,RIR);
      /* calculate ERC luminosity */
      printf("Calculating LERCIR ... ");
      fflush(stdout);
      /* allocate memory */
      vobs = dvector(1,NVERC);
      Lvobs = dvector(1,NVERC);
      theta = (THETAOBS<THETAJ) ? 0.0 : (THETAOBS-THETAJ);
      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
      vobsMAX = 10.0*DSQR(Doppler*gamma[nGAMMA])*epERCIR[NEERC]*mec2h;
      if (vobsMAX>100.0*gamma[nGAMMA]*mec2h)
	vobsMAX = 100.0*gamma[nGAMMA]*mec2h;
      theta = THETAOBS+THETAJ;
      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
      vobsMIN = DSQR(Doppler*gamma[1])*epERCIR[1]*mec2h;
      for(i=1;i<=NVERC;i++)
	vobs[i]=vobsMIN*pow(vobsMAX/vobsMIN,(i-1.0)/(NVERC-1.0));
      N = dvector(1,nGAMMA);
      dtheta = 2.0*THETAJ/NMI;
      for(k=0;k<=NTIME;k++){
	for(i=1;i<=NVERC;i++){
	  Lvobs[i] = 0.0;
	  theta = THETAOBS+THETAJ;
	  do{
	    robs = findr(tobs[k],theta);
	    if(robs>r[1] && robs<r[nRAD]){
	      Doppler = 1.0/(GAMMA*(1.0-cos(theta)*beta(GAMMA)));
	      vpobs = vobs[i]/Doppler;
	      p = locate(r,nRAD,robs);
	      // interpolate Ng into N
	      for(l=1;l<=nGAMMA;l++)
		if (Ng[p][l]>0.0 && Ng[p+1][l]>0.0)
		  N[l] = exp(linear(log(robs),log(r[p]),log(Ng[p][l]),log(r[p+1]),log(Ng[p+1][l])));
		else
		  N[l] = 0.0;
	      ercIR(robs,epERCIR,upeERCIR,NEERC,GAMMA*VIR,4.0/3.0*GAMMA*GAMMA*UIR,RIR);
	      Lvp = LERCIRv(vpobs,gamma,N,nGAMMA,theta);
	      if(theta<(THETAJ-THETAOBS) || THETAOBS==0.0)
		dfi = M_PI;
	      else {
		if (fabs(((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS))))>1.0e0)
		  dfi = 2.0e0*M_PI;
		else
		  dfi = 2.0e0*acos((cos(THETAJ)-cos(theta)*cos(THETAOBS))/(sin(theta)*sin(THETAOBS)));
	      }
	      Lvobs[i] += dfi*Lvp*pow(Doppler,3)*fabs(sin(theta));
	    }
	    theta -= dtheta;
	  } while(theta >= (THETAOBS-THETAJ));
	  Lvobs[i] *= dtheta/(2.0*M_PI*(1.0-cos(THETAJ)));
	}
	save(vobs,Lvobs,"LERCIR",tobs[k]/tobs[0],NVERC);
      }
      /* free memory */
      free_dvector(N,1,nGAMMA);
      free_dvector(Lvobs,1,NVERC);
      free_dvector(vobs,1,NVERC);
      free_dvector(upeERC,1,NEERC);
      free_dvector(epERC,1,NEERC);
      printf("done\n");
    }
    /* free memory */
    free_dmatrix(Ng,1,nRAD,1,nGAMMA);
    free_dvector(gamma,1,nGAMMA);
  }
  free_dvector(r,1,nRAD);


#define NV 500
  double *vS,*LvS,*vSSC,*LvSSC,*vERC,*LvERC,*vERCIR,*LvERCIR;
  char *tf;
  double *v,*L;
  double vmin,vmax;
  extern int locate(double xx[],unsigned long n,double x);
  extern double linear(double x,double x1,double y1,double x2,double y2);
  
  printf("\nAdding spectra ... ");
  filename = (char *)malloc((size_t)(20*sizeof(char)));
  tf = (char *)malloc((size_t)(20*sizeof(char)));
  v = dvector(1,NV);
  L = dvector(1,NV);
  if(SYN){
    vS = dvector(1,NVS);
    LvS = dvector(1,NVS);
  }
  if(SSC){
    vSSC = dvector(1,NVSSC);
    LvSSC = dvector(1,NVSSC);
  }
  if(ERC){
    vERC = dvector(1,NVERC);
    LvERC = dvector(1,NVERC);
  }
  if(ERCIR){
    vERCIR = dvector(1,NVERC);
    LvERCIR = dvector(1,NVERC);
  }
  for(i=0;i<=NTIME;i++){
    vmin = 1.0e20;
    vmax = 0.0;
    sprintf(tf,"%.6f",tobs[i]/tobs[0]);
    if(SYN){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LS_");
      strcat(filename,tf);
      input = fopen(filename,"r");
      for(p=1;p<=NVS;p++) {
        fscanf(input,"%f %f",&tempv,&tempL);
        vS[p] = pow(10.0,(double)tempv);
        LvS[p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
      if(vS[1]<vmin)
	vmin = vS[1];
      if(vS[NVS]>vmax)
	vmax = vS[NVS];
    }
    if(SSC){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LSSC_");
      strcat(filename,tf);
      input = fopen(filename,"r");
      for(p=1;p<=NVSSC;p++) {
        fscanf(input,"%f %f",&tempv,&tempL);
        vSSC[p] = pow(10.0,(double)tempv);
        LvSSC[p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
      if(vSSC[1]<vmin)
	vmin = vSSC[1];
      if(vSSC[NVSSC]>vmax)
	vmax = vSSC[NVSSC];
    }
    if(ERC){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LERC_");
      strcat(filename,tf);
      input = fopen(filename,"r");
      for(p=1;p<=NVERC;p++) {
        fscanf(input,"%f %f",&tempv,&tempL);
        vERC[p] = pow(10.0,(double)tempv);
        LvERC[p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
      if(vERC[1]<vmin)
	vmin = vERC[1];
      if(vERC[NVERC]>vmax)
	vmax = vERC[NVERC];
    }
    if(ERCIR){
      strcpy(filename,"./");
      strcat(filename,INPUT_FILE_ID);
      strcat(filename,"/");
      strcat(filename,"LERCIR_");
      strcat(filename,tf);
      input = fopen(filename,"r");
      for(p=1;p<=NVERC;p++) {
        fscanf(input,"%f %f",&tempv,&tempL);
        vERCIR[p] = pow(10.0,(double)tempv);
        LvERCIR[p] = (tempL==0.0) ? 0.0 : pow(10.0,(double)tempL);
      }
      fclose(input);
      if(vERCIR[1]<vmin)
	vmin = vERCIR[1];
      if(vERCIR[NVERC]>vmax)
	vmax = vERCIR[NVERC];
    }
    /* summation */
    for(k=1;k<=NV;k++){
      v[k] = vmin*pow(vmax/vmin,(k-1.0)/(NV-1.0));
      L[k] = 0.0;
      if(SYN && v[k]>vS[1] && v[k]<vS[NVS]){
	p = locate(vS,NVS,v[k]);
	L[k] += linear(v[k],vS[p],LvS[p],vS[p+1],LvS[p+1]);
      }
      if(SSC && v[k]>vSSC[1] && v[k]<vSSC[NVSSC]){
	p = locate(vSSC,NVSSC,v[k]);
	L[k] += linear(v[k],vSSC[p],LvSSC[p],vSSC[p+1],LvSSC[p+1]);
      }
      if(ERC && v[k]>vERC[1] && v[k]<vERC[NVERC]){
	p = locate(vERC,NVERC,v[k]);
	L[k] += linear(v[k],vERC[p],LvERC[p],vERC[p+1],LvERC[p+1]);
      }
      if(ERCIR && v[k]>vERCIR[1] && v[k]<vERCIR[NVERC]){
	p = locate(vERCIR,NVERC,v[k]);
	L[k] += linear(v[k],vERCIR[p],LvERCIR[p],vERCIR[p+1],LvERCIR[p+1]);
      }
    }
    save(v,L,"L",tobs[i]/tobs[0],NV);
  }
  /* free memory */
  free_dvector(L,1,NV);
  free_dvector(v,1,NV);
  if(ERCIR){
    free_dvector(LvERCIR,1,NVERC);
    free_dvector(vERCIR,1,NVERC);
  }
  if(ERC){
    free_dvector(LvERC,1,NVERC);
    free_dvector(vERC,1,NVERC);
  }
  if(SSC){
    free_dvector(LvSSC,1,NVSSC);
    free_dvector(vSSC,1,NVSSC);
  }
  if(SYN){
    free_dvector(LvS,1,NVS);
    free_dvector(vS,1,NVS);
  }
  printf("done.\n");

  free(tf);
  free(filename);
  free_dvector(tobs,0,NTIME);
  return 0;
}

double
findr(double t, double theta)
     /* finds the radius for a given observed time t>0 - const GAMMA*/
{
  return(LIGHT_SPEED*beta(GAMMA)*t/(1.0-beta(GAMMA)*cos(theta)));
}

/* external Compton luminosity */
#define constC1 (1.21e-34)  // 3/4 sigTch/mc2

double
LERCv(double const v, double const g[], double const N[], int const n,
      double const theta)
{
  double dg,de,Lv;
  double Int1;
  int i,k;
  double f(double g, double e0, double e, double miu);
  double const a1=3.0/8.0,a2=7.0/6.0,a3=23.0/24.0;
  double e,miu;

  e = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  miu = -(cos(theta)-beta(GAMMA))/(1.0-beta(GAMMA)*cos(theta));

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
  Lv  = a1*upeERC[1]/epERC[1]*Int1*dg;

  Int1  = a1*f(g[1],epERC[2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[2],e,miu)*N[n]/g[n];
  Lv  += a2*upeERC[2]/epERC[2]*Int1*dg;

  Int1  = a1*f(g[1],epERC[3],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[3],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[3],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[3],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[3],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[3],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[3],e,miu)*N[n]/g[n];
  Lv  += a3*upeERC[3]/epERC[3]*Int1*dg;

  for(i=4;i<NEERC-2;i++) {
    Int1  = a1*f(g[1],epERC[i],e,miu)*N[1]/g[1];
    Int1 += a2*f(g[2],epERC[i],e,miu)*N[2]/g[2];
    Int1 += a3*f(g[3],epERC[i],e,miu)*N[3]/g[3];
    for(k=4;k<n-2;k++)
      Int1 += f(g[k],epERC[i],e,miu)*N[k]/g[k];
    Int1 += a3*f(g[n-2],epERC[i],e,miu)*N[n-2]/g[n-2];
    Int1 += a2*f(g[n-1],epERC[i],e,miu)*N[n-1]/g[n-1];
    Int1 += a1*f(g[n],epERC[i],e,miu)*N[n]/g[n];
    Lv += upeERC[i]/epERC[i]*Int1*dg;
  }

  Int1  = a1*f(g[1],epERC[NEERC-2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC-2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC-2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC-2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC-2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC-2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC-2],e,miu)*N[n]/g[n];
  Lv  += a3*upeERC[NEERC-2]/epERC[NEERC-2]*Int1*dg;

  Int1  = a1*f(g[1],epERC[NEERC-1],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC-1],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC-1],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC-1],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC-1],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC-1],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC-1],e,miu)*N[n]/g[n];
  Lv  += a2*upeERC[NEERC-1]/epERC[NEERC-1]*Int1*dg;

  Int1  = a1*f(g[1],epERC[NEERC],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERC[NEERC],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERC[NEERC],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERC[NEERC],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERC[NEERC],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERC[NEERC],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERC[NEERC],e,miu)*N[n]/g[n];
  Lv  += a1*upeERC[NEERC]/epERC[NEERC]*Int1*dg;
  
  Lv *= constC1*e*de;
  return Lv;
}

double
LERCIRv(double const v, double const g[], double const N[], int const n,
      double const theta)
{
  double dg,de,Lv;
  double Int1;
  int i,k;
  double f(double g, double e0, double e, double miu);
  double const a1=3.0/8.0,a2=7.0/6.0,a3=23.0/24.0;
  double e,miu;

  e = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  miu = -(cos(theta)-beta(GAMMA))/(1.0-beta(GAMMA)*cos(theta));

  dg = log(g[2]/g[1]);
  de = log(epERCIR[2]/epERCIR[1]);
  Int1  = a1*f(g[1],epERCIR[1],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[1],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[1],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[1],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[1],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[1],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[1],e,miu)*N[n]/g[n];
  Lv  = a1*upeERCIR[1]/epERCIR[1]*Int1*dg;

  Int1  = a1*f(g[1],epERCIR[2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[2],e,miu)*N[n]/g[n];
  Lv  += a2*upeERCIR[2]/epERCIR[2]*Int1*dg;

  Int1  = a1*f(g[1],epERCIR[3],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[3],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[3],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[3],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[3],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[3],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[3],e,miu)*N[n]/g[n];
  Lv  += a3*upeERCIR[3]/epERCIR[3]*Int1*dg;

  for(i=4;i<NEERC-2;i++) {
    Int1  = a1*f(g[1],epERCIR[i],e,miu)*N[1]/g[1];
    Int1 += a2*f(g[2],epERCIR[i],e,miu)*N[2]/g[2];
    Int1 += a3*f(g[3],epERCIR[i],e,miu)*N[3]/g[3];
    for(k=4;k<n-2;k++)
      Int1 += f(g[k],epERCIR[i],e,miu)*N[k]/g[k];
    Int1 += a3*f(g[n-2],epERCIR[i],e,miu)*N[n-2]/g[n-2];
    Int1 += a2*f(g[n-1],epERCIR[i],e,miu)*N[n-1]/g[n-1];
    Int1 += a1*f(g[n],epERCIR[i],e,miu)*N[n]/g[n];
    Lv += upeERCIR[i]/epERCIR[i]*Int1*dg;
  }

  Int1  = a1*f(g[1],epERCIR[NEERC-2],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[NEERC-2],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[NEERC-2],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[NEERC-2],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[NEERC-2],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[NEERC-2],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[NEERC-2],e,miu)*N[n]/g[n];
  Lv  += a3*upeERCIR[NEERC-2]/epERCIR[NEERC-2]*Int1*dg;

  Int1  = a1*f(g[1],epERCIR[NEERC-1],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[NEERC-1],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[NEERC-1],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[NEERC-1],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[NEERC-1],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[NEERC-1],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[NEERC-1],e,miu)*N[n]/g[n];
  Lv  += a2*upeERCIR[NEERC-1]/epERCIR[NEERC-1]*Int1*dg;

  Int1  = a1*f(g[1],epERCIR[NEERC],e,miu)*N[1]/g[1];
  Int1 += a2*f(g[2],epERCIR[NEERC],e,miu)*N[2]/g[2];
  Int1 += a3*f(g[3],epERCIR[NEERC],e,miu)*N[3]/g[3];
  for(k=4;k<n-2;k++)
    Int1 += f(g[k],epERCIR[NEERC],e,miu)*N[k]/g[k];
  Int1 += a3*f(g[n-2],epERCIR[NEERC],e,miu)*N[n-2]/g[n-2];
  Int1 += a2*f(g[n-1],epERCIR[NEERC],e,miu)*N[n-1]/g[n-1];
  Int1 += a1*f(g[n],epERCIR[NEERC],e,miu)*N[n]/g[n];
  Lv  += a1*upeERCIR[NEERC]/epERCIR[NEERC]*Int1*dg;
  
  Lv *= constC1*e*de;
  return Lv;
}

#undef contC1

double
f(double g, double e0, double e, double miu)
{
  double w,wp,b,t,fx;
  
  w = e/g;
  b = 2.0*(1.0-miu)*e0*g;
  if (e>e0 && e<b*g/(1+b)) {
    wp = 1.0-w;
    fx = (1.0+0.5*w*w/wp-2.0*w/(b*wp)+2.0*w*w/(b*b*wp*wp));
    if (fx>0.0)
      return fx;
    else
      return 0.0;
  } else
    return 0.0;
}
