/** @file model.h
    @brief Description of the parameters
    
    This file contains the declarations of parameters and declaration of the function to read them from parameters.dat file.	
*/

#ifndef _MODEL_H
#define _MODEL_H 1

/* energy loses switches */
int ADI,SYN,SSC,ERC,ERCIR,KN;

/* synchrotron self absorption switch */
int SABS;

/* adiabatic cooling constant */
double ADIABG;

/* geometry parameters */
double R0,RINMAX,RMAX,GEOMK,THETAJ,THETAOBS;

/* electrons parameters */
double ELEP1,ELEP2,gammaMIN,gammaM,gammaMAX;

/* time parameters */
double TAVMIN,TAVMAX;

/* bulk Lorentz factor */
double GAMMA;
#define beta(x) (sqrt(1.0-1.0/(x*x)))

/* injected electrons normalization */
double ELEK;

/* magnetic field */
double B0,B1;

/* external radiation */
double RBEL,UBEL,VBEL,BELALPHA,RIR,UIR,VIR;

/* frequency of observation */
double VOBS;

/*injection disspersion*/
double RSIGMA;

/*ERC and B index*/
double KERC;
double KB;

/* maximum injection switch */
double RM;

/* gauss inj switch */
int QGAUSS;

/* triangle inj switch */
int QTRI;

/*input files name suffix*/
char *INPUT_FILE_ID;

/* reading parameters procedure */
int input_parameters(char *id);

/* magnetic field */
double magB( double x);

#endif /* _MODEL_H */
