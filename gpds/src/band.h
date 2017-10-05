#ifndef BAND_H
#define BAND_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>

using namespace std;



extern "C" {

  void bmatvecmult(double *a, double *b, int *bandsize, int *matdim, double *result);

  void bmatvecmultT(double *a, double *b, int *bandsize, int *matdim, double *result);


  void xthetallikBandC( double *xtheta, double *Vmphi, double *VKinv, double *VCinv,
                        double *Rmphi, double *RKinv, double *RCinv, int *bandsize, int *nn,
                        double *sigma, double *yobs, double *ret, double *retgrad);

}

// g++ band.cpp -o band.o -lopenblas -llapack -lm -Wall -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include
#endif
