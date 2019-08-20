#include "mex.h"
#include "armaMex.hpp"
#include "../gpds_cpp/tgtdistr.h"
#include "../gpds_cpp/dynamicalSystemModels.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//     Rcpp::List xthetallikRcpp(const arma::mat & yobs, 
//                           const Rcpp::List & covAllDimInput,
//                           const arma::vec & sigmaInput,
//                           const arma::vec & initial,
//                           const std::string modelName = "FN",
//                           const bool useBand = false,
//                           const Rcpp::NumericVector & priorTemperatureInput = 1.0){
  
   mat yobs = armaGetPr(prhs[0]);
   //int nfields = mxGetNumberOfFields(prhs[1]);
   int NStructElems = mxGetNumberOfElements(prhs[1]);  // need this one
   
   //mat Ctest = armaGetPr(mxGetField(prhs[1], 0, "C"));
   //Ctest.print();
   vec sigmaInput = armaGetPrVec(prhs[2]);
   vec initial = armaGetPrVec(prhs[3]);

   //cout << nfields << " " << NStructElems;
   
   char *str1;
   str1 = mxArrayToString(prhs[4]);
   string modelName = str1;
   
   
   bool useBand = mxGetScalar(prhs[5]);
   
   arma::vec priorTemperature = {1.0};
    
  //const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  vec sigma( yobs.n_cols);
  if(sigmaInput.size() == 1){
    sigma.fill( as_scalar( sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput size not right");
  }
  
  
  
  vector<gpcov> covAllDimensions(yobs.n_cols);
  for(unsigned j = 0; j < yobs.n_cols; j++){
    //covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
    covAllDimensions[j].Cinv = armaGetPr(mxGetField(prhs[1], j, "Cinv"));
    covAllDimensions[j].mphi = armaGetPr(mxGetField(prhs[1], j, "mphi"));
    covAllDimensions[j].Kinv = armaGetPr(mxGetField(prhs[1], j, "Kinv"));
    covAllDimensions[j].CinvBand = armaGetPr(mxGetField(prhs[1], j, "CinvBand"));
    covAllDimensions[j].mphiBand = armaGetPr(mxGetField(prhs[1], j, "mphiBand"));
    covAllDimensions[j].KinvBand = armaGetPr(mxGetField(prhs[1], j, "KinvBand"));
    covAllDimensions[j].mu = armaGetPrVec(mxGetField(prhs[1], j, "mu"));
    covAllDimensions[j].dotmu = armaGetPrVec(mxGetField(prhs[1], j, "dotmu"));

  }
  
//   const Rcpp::NumericMatrix & Cinv = as<const NumericMatrix>(cov_r["Cinv"]);
//   const Rcpp::NumericMatrix & mphi = as<const NumericMatrix>(cov_r["mphi"]);
//   const Rcpp::NumericMatrix & Kinv = as<const NumericMatrix>(cov_r["Kinv"]);
//   const Rcpp::NumericMatrix & CinvBand = as<const NumericMatrix>(cov_r["CinvBand"]);
//   const Rcpp::NumericMatrix & mphiBand = as<const NumericMatrix>(cov_r["mphiBand"]);
//   const Rcpp::NumericMatrix & KinvBand = as<const NumericMatrix>(cov_r["KinvBand"]);
//   const Rcpp::NumericVector & mu = as<const NumericVector>(cov_r["mu"]);
//   const Rcpp::NumericVector & dotmu = as<const NumericVector>(cov_r["dotmu"]);  
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "Hes1-log"){
    model = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  lp ret = xthetallik(initial, covAllDimensions, sigma, yobs, model, useBand, priorTemperature);
  
  plhs[0] = mxCreateDoubleScalar(ret.value);
  plhs[1] = mxCreateDoubleMatrix(ret.gradient.n_rows, ret.gradient.n_cols, mxREAL); 
  armaSetPr(plhs[1], ret.gradient);
 
  
  
  //cout << ret.value << "\n" << ret.gradient;
  //return List::create(Named("value")=ret.value,
  //                    Named("grad")=ret.gradient);
    
}