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
  const mat thetaInitExogenous= armaGetPr(prhs[6]);
  const mat muExogenous= armaGetPr(prhs[7]);
  const mat dotmuExogenous= armaGetPr(prhs[8]);

  const double priorTemperatureLevel = mxGetScalar(prhs[9]);
  const double priorTemperatureDeriv = mxGetScalar(prhs[10]);
  const double priorTemperatureObs = mxGetScalar(prhs[11]);

  char *str1;
  str1 = mxArrayToString(prhs[12]);
  string kernel = str1;
  
  const int nstepsHmc = mxGetScalar(prhs[13]);
  const double burninRatioHmc = mxGetScalar(prhs[14]);
  const unsigned int niterHmc = mxGetScalar(prhs[15]);
  const double stepSizeFactorHmc = mxGetScalar(prhs[16]);
  const int nEpoch = mxGetScalar(prhs[17]);
  const int bandSize = mxGetScalar(prhs[18]);
  bool useFrequencyBasedPrior = mxGetLogicals(prhs[19])[0];
  bool useBand = mxGetLogicals(prhs[20])[0];
  bool useMean = mxGetLogicals(prhs[21])[0];
  bool useScalerSigma = mxGetLogicals(prhs[22])[0];
  bool useFixedSigma = mxGetLogicals(prhs[23])[0];
  bool verbose = mxGetLogicals(prhs[24])[0];    

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
      mxDestroyArray(x_matlab);
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
          thetaInitExogenous,
          muExogenous,
          dotmuExogenous,
          priorTemperatureLevel,
          priorTemperatureDeriv,
          priorTemperatureObs,
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
  
  solver.initMissingComponent();
  solver.sampleInEpochs();
  
  cube ret = solver.llikxthetasigmaSamples;
  mat retphi = solver.phiAllDimensions;

  mwSize dimensions[3];
  dimensions[0] = ret.n_rows;
  dimensions[1] = ret.n_cols;
  dimensions[2] = ret.n_slices;
  
  mwSize phidim[2];
  phidim[0] = retphi.n_rows;
  phidim[1] = retphi.n_cols;

  plhs[0] = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2, phidim, mxDOUBLE_CLASS, mxREAL);
  armaSetCubePr(plhs[0], ret);
  armaSetPr(plhs[1], retphi);
  

}
