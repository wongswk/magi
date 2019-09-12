#include "mex.h"
#include "armaMex.hpp"
#include "../../gpds_cpp/GpdsSolver.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
  const mat yFull = armaGetPr(prhs[0]);
  
  const vec tvecFull = armaGetPrVec(prhs[2]);
  const vec sigmaExogenous = armaGetPrVec(prhs[3]);
  const mat phiExogenous = armaGetPr(prhs[4]);
  const mat xInitExogenous= armaGetPr(prhs[5]);
  const mat muExogenous= armaGetPr(prhs[6]);
  const mat dotmuExogenous= armaGetPr(prhs[7]);

  const double priorTemperatureLevel = mxGetScalar(prhs[8]);
  const double priorTemperatureDeriv = mxGetScalar(prhs[9]);

  char *str1;
  str1 = mxArrayToString(prhs[10]);
  string kernel = str1;
  
  const int nstepsHmc = mxGetScalar(prhs[11]);
  const double burninRatioHmc = mxGetScalar(prhs[12]);
  const unsigned int niterHmc = mxGetScalar(prhs[13]);
  const double stepSizeFactorHmc = mxGetScalar(prhs[14]);
  const int nEpoch = mxGetScalar(prhs[15]);
  const int bandSize = mxGetScalar(prhs[16]);
  bool useFrequencyBasedPrior = mxGetLogicals(prhs[17])[0];
  bool useBand = mxGetLogicals(prhs[18])[0];
  bool useMean = mxGetLogicals(prhs[19])[0];
  bool useScalerSigma = mxGetLogicals(prhs[20])[0];
  bool useFixedSigma = mxGetLogicals(prhs[21])[0];
  bool verbose = mxGetLogicals(prhs[22])[0];    

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
      mat temp_ret = armaGetPr(lhs_f[0]);
      mat ret = temp_ret;
    
      mxDestroyArray(x_matlab);
      mxDestroyArray(theta_matlab);
      mxDestroyArray(lhs_f[0]);      
      
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
      cube temp_ret = armaGetCubePr(lhs_f[0]);

      cube ret = temp_ret;
      mxDestroyArray(x_matlab);
      mxDestroyArray(theta_matlab);
      mxDestroyArray(lhs_f[0]);

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
      cube temp_ret = armaGetCubePr(lhs_f[0]);

      cube ret = temp_ret;
      mxDestroyArray(theta_matlab);
      mxDestroyArray(lhs_f[0]);           
      
      return ret;
  };
  
  OdeSystem model;
  model = OdeSystem(fOde, fOdeDx, fOdeDtheta, lb, ub);

  GpdsSolver solver(yFull,
          model,
          tvecFull,
          sigmaExogenous,
          phiExogenous,
          xInitExogenous,
          muExogenous,
          dotmuExogenous,
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
  
  solver.sampleInEpochs();
  
  mat ret = join_vert(solver.llikSamples.t(), solver.xthetasigmaSamples);
  
  plhs[0] = mxCreateDoubleMatrix(ret.n_rows, ret.n_cols, mxREAL); 
  armaSetPr(plhs[0], ret);
  

}
