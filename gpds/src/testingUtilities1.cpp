#include<iostream>
#include "band.h"

// [[Rcpp::export]]
int bandTest(){
  std::cout << "in bandTest\n";
  return mainBand();
}
