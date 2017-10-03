#ifndef BANDMATLIK_H
#define BANDMATLIK_H

#include <R.h>
#include <R_ext/BLAS.h>
#include <stdio.h>
#include <iostream>

using namespace std;


extern "C" {
  
  void bmatvecmult(double *a, double *b, int *bandsize, int *matdim, double *result);
  
  void bmatvecmultT(double *a, double *b, int *bandsize, int *matdim, double *result);   
  
  void xthetallik( double *xtheta, double *Vmphi, double *VKinv, double *VCinv,
                   double *Rmphi, double *RKinv, double *RCinv, int *bandsize, int *nn,
                   double *sigma, double *yobs, double *ret, double *retgrad);
  
}
#endif