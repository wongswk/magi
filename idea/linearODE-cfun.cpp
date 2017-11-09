// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
double logPosteriorC( const vec & ab,
                      const mat & AxInvPartConst,
                      const mat & Kinv,
                      const mat & AxInvPartConstKinvMphiSum,
                      const vec & muxRawPartConst,
                      const vec & KinvRowSum,
                      const vec & yContribution,
                      const double & KinvSum) {
  
  mat AxInv =  AxInvPartConst + std::pow(ab(1), 2) * Kinv - ab(1) * AxInvPartConstKinvMphiSum;
  
  vec muxRaw = ab(0) * muxRawPartConst - ab(0) * ab(1) * KinvRowSum + yContribution;
  vec mux = solve( AxInv, muxRaw);
  
  double llik;
  double sign;
  
  log_det(llik, sign, AxInv); 
  
  
  llik -= std::pow(ab(0), 2) * KinvSum - accu( muxRaw % mux);
  llik *= 0.5;
  return llik;
}

