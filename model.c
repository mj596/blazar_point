#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include <string.h>
#include <sys/stat.h>

int
input_parameters(char* id)
{
  FILE *input;
  int i,tmpi;
  float tmpx;
  double tmpy;
  char tmps[256];
  char *filename;

  filename = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(filename,"parameters");
  strcat(filename,"_");
  strcat(filename,id);
  strcat(filename,".dat");

  INPUT_FILE_ID = (char *)malloc((size_t)(20*sizeof(char)));
  strcpy(INPUT_FILE_ID,id);

  /*make a data directory*/
  mkdir(id,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  printf("\nReading parameters from... %s\n",filename);

  if((input=fopen(filename,"r"))==NULL){
    fprintf(stderr,"error: can't open %s\n",filename);
    exit(1);
  }
  
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  ADI = (int)tmpi;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  SYN = (int)tmpi;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  SSC = (int)tmpi;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  ERC = (int)tmpi;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  ERCIR = (int)tmpi;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  ADIABG = (double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  SABS = (int)tmpx;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  KN = (int)tmpi;
  printf("Cooling switches:\n\tADI %d  SYN %d  SSC %d  ERC %d  ERCIR %d\n\tADIABG=%4.2f  SABS %d  KN %d\n",ADI,SYN,SSC,ERC,ERCIR,ADIABG,SABS,KN);
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  gammaMIN=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  gammaM=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  gammaMAX=(double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  ELEP1 = (double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  ELEP2 = (double)tmpx;
  fscanf(input,"%lf%[^\n]\n",&tmpy,tmps);
  ELEK = (double)tmpy;
  printf("Electrons:\n\tgammaMIN=%5.1e  gammaM=%5.1e gammaMAX=%5.1e\n\tELEP1=%3.1f ELEP2=%3.1f ELEK=%6.2e\n",gammaMIN,gammaM,gammaMAX,ELEP1,ELEP2,ELEK);
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  GAMMA=(double)tmpx;
  printf("Bulk Lorentz factor:\n\tGAMMA=%5.1f\n",GAMMA);
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  THETAJ=(double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  THETAOBS=(double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  GEOMK=(double)tmpx;
  printf("Geometry:\n\tTHETAJ=%4.2f  THETAOBS=%4.2f  GEOMK=%3.1f\n",THETAJ,THETAOBS,GEOMK);
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  R0=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  RINMAX = (double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  RMAX=(double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  TAVMIN = (double)tmpx;
  fscanf(input,"%f%[^\n]\n",&tmpx,tmps);
  TAVMAX = (double)tmpx;
  printf("Observations:\n\tR0=%6.1e  RINMAX=%6.1e RMAX=%6.1e\n\tTAVMIN=%6.1e  TAVMAX=%6.1e\n",R0,RINMAX,RMAX,TAVMIN,TAVMAX);
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  B0=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  B1=(double)tmpx;
  printf("Magnetic field:\n\tB0=%6.1e  B1=%6.1e\n",B0,B1);
  fscanf(input,"%le%[^\n]\n",&tmpy,tmps);
  RBEL = (double)tmpy;
  fscanf(input,"%le%[^\n]\n",&tmpy,tmps);
  UBEL =(double)tmpy;
  fscanf(input,"%le%[^\n]\n",&tmpy,tmps);
  VBEL = (double)tmpy;
  fscanf(input,"%le%[^\n]\n",&tmpy,tmps);
  BELALPHA = (double)tmpy;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  RIR = (double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  UIR = (double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  VIR = (double)tmpx;
  printf("External fields:\n\tVBEL=%6.1e  RBEL=%6.1e  UBEL=%6.1e  BELALPHA=%6.1e\n\tVIR=%6.1e  RIR=%6.1e  UIR=%6.1e\n",VBEL,RBEL,UBEL,BELALPHA,VIR,RIR,UIR);
  printf("Observing frequency:\n\t");
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  VOBS = (double)tmpx;
  printf("VOBS = %4.2e  \n",VOBS);
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  RSIGMA=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  KERC=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  KB=(double)tmpx;
  fscanf(input,"%e%[^\n]\n",&tmpx,tmps);
  RM = (double)tmpx;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  QGAUSS = (int)tmpi;
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  QTRI = (int)tmpi;
  printf("Extra:\n\tRSIGMA=%4.2e  KERC=%2.1f  KB=%2.1f RM=%4.2e  QGAUSS=%d QTRI=%d\n",RSIGMA,KERC,KB,RM,QGAUSS,QTRI);
  fscanf(input,"%d%[^\n]\n",&tmpi,tmps);
  printf("\n");
  printf("... done.\n");
  return 0;
}

double magB( double x )
{
  return B0+(RM*B1)/x;
}
