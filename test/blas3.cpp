#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

int main( int argc, char** argv ){
  // you can define the arrays in one of two ways
  // on the heap
  double *a = (double*) malloc( 3 * sizeof(double) );
  a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
  // on the stack
  double b[3] = { 4.0, 5.0, 6.0 };

  double dot_product = cblas_ddot( 3, a, 1, b, 1 );
  printf(" The dot product is: %f \n",dot_product );

  return 0;
};
//  g++ blas3.cpp -o lap.o -lopenblas -llapack -lm -Wall -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include