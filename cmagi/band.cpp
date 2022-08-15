#include "band.h"

typedef size_t blasint;

extern "C" {

void dgbmv_(const char *, const blasint *, const blasint *, const blasint *, const blasint *, const double *, const double *, const blasint *,
            const double *, const blasint *, const double *, double *, const blasint *);

void bmatvecmult(const double *a, const double *b, const int *bandsize, const int *matdim, double *result) {

    double zero = 0.0;
    double one = 1.0;
    blasint ione = 1;
    char nthechar = 'n';
    //std::cout << *a << " " << *bandsize << " " << *matdim <<  endl;
    blasint nrow = 2 * *bandsize +  1;
    blasint bmatdim = *matdim;
    blasint bbandsize = *bandsize;
    //std::cout << nrow << endl;

    dgbmv_(&nthechar, &bmatdim, &bmatdim, &bbandsize, &bbandsize, &one, a, &nrow, b, &ione, &zero, result,
           &ione);


}

void bmatvecmultT(const double *a, const double *b, const int *bandsize, const int *matdim, double *result) {

    double zero = 0.0;
    double one = 1.0;
    blasint ione = 1;
    blasint nrow = 2 * *bandsize +  1;
    char tthechar = 't';
    blasint bmatdim = *matdim;
    blasint bbandsize = *bandsize;

    dgbmv_(&tthechar,  &bmatdim, &bmatdim, &bbandsize, &bbandsize, &one, a, &nrow, b, &ione, &zero, result,
           &ione);

}

}



// g++ band.cpp -o band.o -lopenblas -llapack -lm -Wall -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include
