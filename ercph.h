/** @file ercph.h
    @brief Energy density of external radiation fields
    
    This file contains the declarations of functions to calculate the energy density of the external radiation fields.	
*/

#ifndef _ERCPH_H
#define _ERCPH_H 1

void
ercPL(double* ep, double* upe, int n,
      double epext, double upext, double alpha);

void
ercBB(double* ep, double* upe, int n,
      double epext, double upext);

void
ercBEL(double r, double* ep, double* upe, int n,
       double epext, double upext, double rbel, double alpha);

void
ercIR(double r, double* ep, double* upe, int n,
      double epext, double upext, double rir);

#endif /* _ERCPH_H */
