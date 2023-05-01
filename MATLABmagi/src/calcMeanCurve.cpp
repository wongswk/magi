#include "mex.h"
#include "armaMex.hpp"
#include "../../cmagi/gpsmoothing.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  const vec tvec = armaGetPrVec(prhs[0]);
  const vec yobs = armaGetPrVec(prhs[1]);
  const vec tnew = armaGetPrVec(prhs[2]);
  const mat phi = armaGetPr(prhs[3]);
  const vec sigma = armaGetPrVec(prhs[4]);

  char *str1;
  str1 = mxArrayToString(prhs[5]);
  string kerneltype = str1;

  bool deriv = mxGetLogicals(prhs[6])[0];
  
  cube ret = calcMeanCurve(tvec, yobs, tnew, phi, sigma, kerneltype, deriv);
  
  mwSize dimensions[3];
  dimensions[0] = ret.n_rows;
  dimensions[1] = ret.n_cols;
  dimensions[2] = ret.n_slices;
  
  plhs[0] = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS, mxREAL);
    
  armaSetPr(plhs[0], ret);
    

}