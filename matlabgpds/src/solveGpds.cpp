#include "mex.h"
#include "armaMex.hpp"
#include "../../gpds_cpp/GpdsSolver.h"
#include "../../gpds_cpp/GpdsSolver.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

// arma::mat solveGpds(const arma::mat & yFull,
//                     const OdeSystem & odeModel,
//                     const arma::vec & tvecFull,
//                     const arma::vec & sigmaExogenous = arma::vec(),
//                     const double priorTemperatureLevel = 1,
//                     const double priorTemperatureDeriv = 1,
//                     std::string kernel = "generalMatern",
//                     const int nstepsHmc = 500,
//                     const double burninRatioHmc = 0.5,
//                     const unsigned int niterHmc = 10000,
//                     const double stepSizeFactorHmc = 1,
//                     const int nEpoch = 10,
//                     const int bandSize = 20,
//                     bool useFrequencyBasedPrior = false,
//                     bool useBand = true,
//                     bool useMean = true,
//                     bool useScalerSigma = false,
//                     bool useFixedSigma = false,
//                     bool verbose = false) 
  
  const mat yFull = armaGetPr(prhs[0]);
  
  const vec tvecFull = armaGetPrVec(prhs[2]);
  const vec sigmaExogenous = armaGetPrVec(prhs[3]);
  const double priorTemperatureLevel = mxGetScalar(prhs[4]);
  const double priorTemperatureDeriv = mxGetScalar(prhs[5]);

  char *str1;
  str1 = mxArrayToString(prhs[6]);
  string kernel = str1;
  
  const int nstepsHmc = mxGetScalar(prhs[7]);
  const double burninRatioHmc = mxGetScalar(prhs[8]);
  const unsigned int niterHmc = mxGetScalar(prhs[9]);
  const double stepSizeFactorHmc = mxGetScalar(prhs[10]);
  const int nEpoch = mxGetScalar(prhs[11]);
  const int bandSize = mxGetScalar(prhs[12]);
  bool useFrequencyBasedPrior = mxGetLogicals(prhs[13])[0];
  bool useBand = mxGetLogicals(prhs[14])[0];
  bool useMean = mxGetLogicals(prhs[15])[0];
  bool useScalerSigma = mxGetLogicals(prhs[16])[0];
  bool useFixedSigma = mxGetLogicals(prhs[17])[0];
  bool verbose = mxGetLogicals(prhs[18])[0];    

  mxArray *fOde_matlab = const_cast<mxArray *>(mxGetField(prhs[1], 0, "fOde"));
  mxArray *fOdeDx_matlab = const_cast<mxArray *>(mxGetField(prhs[1], 0, "fOdeDx"));
  mxArray *fOdeDtheta_matlab = const_cast<mxArray *>(mxGetField(prhs[1], 0, "fOdeDtheta"));
  
  const vec lb = armaGetPrVec(mxGetField(prhs[1], 0, "thetaLowerBound"));
  const vec ub = armaGetPrVec(mxGetField(prhs[1], 0, "thetaUpperBound"));
  
  if( !mxIsClass( fOde_matlab, "function_handle")) {
    mexErrMsgTxt("fOde is not a function handle.");
  }

  if( !mxIsClass( fOdeDx_matlab, "function_handle")) {
    mexErrMsgTxt("fOdeDx is not a function handle.");
  }

  if( !mxIsClass( fOdeDtheta_matlab, "function_handle")) {
    mexErrMsgTxt("fOdeDtheta is not a function handle.");
  }
    
  const std::function<mat (vec, mat)> & fOde = [& fOde_matlab](const vec & theta, const mat & x) -> mat {
      mxArray * theta_matlab = mxCreateDoubleMatrix(theta.size(), 1, mxREAL);
      armaSetPr(theta_matlab, theta);
      mxArray * x_matlab = mxCreateDoubleMatrix(x.n_rows, x.n_cols, mxREAL);
      armaSetPr(x_matlab, x);
      
      mxArray *rhs_f[3];
      rhs_f[0] = fOde_matlab;
      rhs_f[1] = theta_matlab;
      rhs_f[2] = x_matlab;
      
      mxArray *lhs_f[1];
      mexCallMATLAB(1,lhs_f,3,rhs_f,"feval");
      mat ret = armaGetPr(lhs_f[0]);
      
      return ret;
  };

  const std::function<cube (vec, mat)> & fOdeDx = [& fOdeDx_matlab](const vec & theta, const mat & x) -> cube {
      mxArray * theta_matlab = mxCreateDoubleMatrix(theta.size(), 1, mxREAL);
      armaSetPr(theta_matlab, theta);
      mxArray * x_matlab = mxCreateDoubleMatrix(x.n_rows, x.n_cols, mxREAL);
      armaSetPr(x_matlab, x);
      
      mxArray *rhs_f[3];
      rhs_f[0] = fOdeDx_matlab;
      rhs_f[1] = theta_matlab;
      rhs_f[2] = x_matlab;
      
      mxArray *lhs_f[1];
      mexCallMATLAB(1,lhs_f,3,rhs_f,"feval");
      cube ret = armaGetCubePr(lhs_f[0]);
      
      return ret;
  };
  
  const std::function<cube (vec, mat)> & fOdeDtheta = [& fOdeDtheta_matlab](const vec & theta, const mat & x) -> cube {
      mxArray * theta_matlab = mxCreateDoubleMatrix(theta.size(), 1, mxREAL);
      armaSetPr(theta_matlab, theta);
      mxArray * x_matlab = mxCreateDoubleMatrix(x.n_rows, x.n_cols, mxREAL);
      armaSetPr(x_matlab, x);
      
      mxArray *rhs_f[3];
      rhs_f[0] = fOdeDtheta_matlab;
      rhs_f[1] = theta_matlab;
      rhs_f[2] = x_matlab;
      
      mxArray *lhs_f[1];
      mexCallMATLAB(1,lhs_f,3,rhs_f,"feval");
      cube ret = armaGetCubePr(lhs_f[0]);
      
      return ret;
  };
  
  OdeSystem model;
  model = OdeSystem(fOde, fOdeDx, fOdeDtheta, lb, ub);

  GpdsSolver         solver(yFull,
                            model,
                            tvecFull,
                            sigmaExogenous,
                            priorTemperatureLevel,
                            priorTemperatureDeriv,
                            kernel,
                            nstepsHmc,
                            burninRatioHmc,
                            niterHmc,
                            stepSizeFactorHmc,
                            nEpoch,
                            bandSize,
                            useFrequencyBasedPrior,
                            useBand,
                            useMean,
                            useScalerSigma,
                            useFixedSigma,
                            verbose);
  
  
  solver.setupPhiSigma();
  if(verbose){
      std::cout << "phi = \n" << solver.phiAllDimensions << "\n";
  }
  solver.initXmudotmu();
  
  solver.initTheta();
  if(verbose){
      std::cout << "thetaInit = \n" << solver.thetaInit << "\n";
  }
  
  solver.doHMC();
  
  mat ret = join_vert(solver.llikSamples.t(), solver.xthetasigmaSamples);
  
  plhs[0] = mxCreateDoubleMatrix(ret.n_rows, ret.n_cols, mxREAL); 
  armaSetPr(plhs[0], ret);
  

}
