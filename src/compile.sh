# need to comment out <RcppArmadillo.h> and namespace Rcpp;
g++ hmc.cpp -o hmc.o -O2 -larmadillo