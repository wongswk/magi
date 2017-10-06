#include "band.h"

typedef int blasint;

extern "C" {
  
  void dgbmv_(const char *, const blasint *, const blasint *, const blasint *, const blasint *, const double *, const double *, const blasint *,
                const double *, const blasint *, const double *, double *, const blasint *);

  void bmatvecmult(const double *a, const double *b, const int *bandsize, const int *matdim, double *result) {

    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    char nthechar = 'n';
    //cout << *a << " " << *bandsize << " " << *matdim <<  endl;
    int nrow = 2 * *bandsize +  1;

    //cout << nrow << endl;

    dgbmv_(&nthechar, matdim, matdim, bandsize, bandsize, &one, a, &nrow, b, &ione, &zero, result,
             &ione);


  }

  void bmatvecmultT(const double *a, const double *b, const int *bandsize, const int *matdim, double *result) {

    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    int nrow = 2 * *bandsize +  1;
    char tthechar = 't';

    dgbmv_(&tthechar, matdim, matdim, bandsize, bandsize, &one, a, &nrow, b, &ione, &zero, result,
             &ione);

  }


  void xthetallikBandC( const double *xtheta, const double *Vmphi, const double *VKinv, const double *VCinv,
                        const double *Rmphi, const double *RKinv, const double *RCinv, const int *bandsize, const int *nn,
                        const double *sigma, const double *yobs, double *ret, double *retgrad) {

    int n = *nn;
    int i,j;

    if (xtheta[2*n] < 0 || xtheta[2*n+1] < 0 || xtheta[2*n+2] < 0) {
      *ret = -1e+9;
      for (i = 0; i < n; i++)
        retgrad[i] = 0.0;
      return;
    }

    double *Vsm, *Rsm, *Vdt, *Rdt, *frV, *frR, *fitLevelErrorV, *fitLevelErrorR;
    double *tempV, *tempR, *tempV2, *tempR2, *tempV3, *tempR3;
    Vsm = new double[n];
    Rsm = new double[n];
    Vdt = new double[n];
    Rdt = new double[n];
    frV = new double[n];
    frR = new double[n];
    fitLevelErrorV = new double[n];
    fitLevelErrorR = new double[n];
    tempV = new double[n];
    tempR = new double[n];
    tempV2 = new double[n];
    tempR2 = new double[n];
    tempV3 = new double[n];
    tempR3 = new double[n];


    double *C1, *VC2, *RC2, *C3;
    C1 = new double[2*n+3];
    VC2 = new double[2*n+3];
    RC2 = new double[2*n+3];
    C3 = new double[2*n+3];


    double theta[3];
    theta[0] = xtheta[2*n];
    theta[1] = xtheta[2*n+1];
    theta[2] = xtheta[2*n+2];

    for (i = 0; i < n; i++) {
      Vsm[i] = xtheta[i];
      Rsm[i] = xtheta[i+n];

      Vdt[i] = theta[2] * (Vsm[i] - pow(Vsm[i],3) / 3.0 + Rsm[i]);
      Rdt[i] = -1.0/theta[2] * ( Vsm[i] - theta[0] + theta[1] * Rsm[i]);
    }

    // V
    bmatvecmult(Vmphi,Vsm,bandsize,nn,frV);

    for (i = 0; i < n; i++) {
      //cout << frV[i] << " ";
      frV[i] = Vdt[i] - frV[i];

      if (!isnan(yobs[i]))
        fitLevelErrorV[i] = Vsm[i] - yobs[i];
      else
        fitLevelErrorV[i] = 0.0;
    }

    double res00 = 0.0;
    for (i = 0; i < n; i++)
      res00 += fitLevelErrorV[i] * fitLevelErrorV[i];
    res00 = -0.5 * res00 / (*sigma * *sigma);

    bmatvecmult(VKinv,frV,bandsize,nn,tempV);
    double res01 = 0.0;
    for (i = 0; i < n; i++)
      res01 += frV[i] * tempV[i];
    res01 *= -0.5;

    bmatvecmult(VCinv,Vsm,bandsize,nn,tempV2);
    double res02 = 0.0;
    for (i = 0; i < n; i++)
      res02 += Vsm[i] * tempV2[i];
    res02 *= -0.5;

    //cout << res00 << " " << res01 << " " << res02 << endl;

    // R
    bmatvecmult(Rmphi,Rsm,bandsize,nn,frR);
    for (i = 0; i < n; i++) {
      //cout << frV[i] << " ";
      frR[i] = Rdt[i] - frR[i];

      if (!isnan(yobs[i+n]))
        fitLevelErrorR[i] = Rsm[i] - yobs[i+n];
      else
        fitLevelErrorR[i] = 0.0;
    }

    double res10 = 0.0;
    for (i = 0; i < n; i++)
      res10 += fitLevelErrorR[i] * fitLevelErrorR[i];
    res10 = -0.5 * res10 / (*sigma * *sigma);

    bmatvecmult(RKinv,frR,bandsize,nn,tempR);
    double res11 = 0.0;
    for (i = 0; i < n; i++)
      res11 += frR[i] * tempR[i];
    res11 *= -0.5;

    bmatvecmult(RCinv,Rsm,bandsize,nn,tempR2);
    double res12 = 0.0;
    for (i = 0; i < n; i++)
      res12 += Rsm[i] * tempR2[i];
    res12 *= -0.5;

    //cout << res10 << " " << res11 << " " << res12 <<  endl;

    // gradient
    // V contrib
    double *Vtemp;
    int m = (2* *bandsize+1)*n;
    Vtemp = new double[m];
    for (i = 0; i < m; i++) {
      Vtemp[i] = -Vmphi[i];
    }

    j = 0;
    for (i = *bandsize; i < m; i+= 2* *bandsize+1) {
      //cout << Vtemp[i] << " ";
      Vtemp[i] += theta[2] * (1 - Vsm[j] * Vsm[j]);
      j++;
    }
    bmatvecmultT(Vtemp,tempV,bandsize,nn,tempV3);
    VC2[2*n] = 0.0;
    VC2[2*n+1] = 0.0;
    VC2[2*n+2] = 0.0;
    for (i=0; i < n; i++) {
      VC2[i] = tempV3[i];
      VC2[i+n] = theta[2] * tempV[i];
      VC2[2*n+2] += tempV[i] * Vdt[i];
    }
    VC2[2*n+2] /= theta[2];

    // R contrib
    double *Rtemp;
    Rtemp = new double[m];
    for (i = 0; i < m; i++) {
      Rtemp[i] = -Rmphi[i];
    }
    for (i = *bandsize; i < m; i+= 2* *bandsize+1) {
      //cout << Vtemp[i] << " ";
      Rtemp[i] -= theta[1]/theta[2];
    }
    bmatvecmultT(Rtemp,tempR,bandsize,nn,tempR3);
    RC2[2*n] = 0.0;
    RC2[2*n+1] = 0.0;
    RC2[2*n+2] = 0.0;
    for (i=0; i < n; i++) {
      RC2[i] =  -tempR[i] / theta[2];
      RC2[i+n] = tempR3[i];

      RC2[2*n] += tempR[i];
      RC2[2*n+1] -= Rsm[i] * tempR[i];
      RC2[2*n+2] -= tempR[i] * Rdt[i];
    }
    RC2[2*n] /= theta[2];
    RC2[2*n+1] /= theta[2];
    RC2[2*n+2] /= theta[2];

    for (i = 0; i < n; i++) {
      C3[i] = tempV2[i];
      C3[i+n] = tempR2[i];
    }
    C3[2*n] = 0.0;
    C3[2*n+1] = 0.0;
    C3[2*n+2] = 0.0;

    for (i = 0; i < n; i++) {
      C1[i] = fitLevelErrorV[i] / ( *sigma * *sigma);
      C1[i+n] = fitLevelErrorR[i] / ( *sigma * *sigma);
    }
    C1[2*n] = 0.0;
    C1[2*n+1] = 0.0;
    C1[2*n+2] = 0.0;

    for (i = 0; i < 2*n+3; i++) {
      retgrad[i] = - (VC2[i] + RC2[i] + C3[i] + C1[i]);
    }

    *ret = res00 + res01 + res02 + res10 + res11 + res12;

    delete[] Vsm;
    delete[] Rsm;
    delete[] Vdt;
    delete[] Rdt;
    delete[] frV;
    delete[] frR;
    delete[] fitLevelErrorV;
    delete[] fitLevelErrorR;
    delete[] tempV;
    delete[] tempR;
    delete[] tempV2;
    delete[] tempR2;
    delete[] tempV3;
    delete[] tempR3;
    delete[] C1;
    delete[] VC2;
    delete[] RC2;
    delete[] C3;
    delete[] Vtemp;
    delete[] Rtemp;

  }

}



// g++ band.cpp -o band.o -lopenblas -llapack -lm -Wall -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include
