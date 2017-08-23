# need to comment out:
#     <RcppArmadillo.h> 
#     namespace Rcpp;
#     R wrapper for basic_hmcC
g++ hmc.cpp -o hmc.o -O2 -larmadillo
./hmc.o