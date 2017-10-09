#ifndef BAND_H
#define BAND_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>

using namespace std;



extern "C" {

  void xthetallikBandC( const double *xtheta, const double *Vmphi, const double *VKinv, const double *VCinv,
                        const double *Rmphi, const double *RKinv, const double *RCinv, const int *bandsize, const int *nn,
                        const double *sigma, const double *yobs, double *ret, double *retgrad);
  // previous .C wrapper doesn't work with Rcpp auto-generated wrapper
}

// g++ band.cpp -o band.o -lopenblas -llapack -lm -Wall -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include
#endif
