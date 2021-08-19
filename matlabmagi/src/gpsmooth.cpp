#include "mex.h"
#include "armaMex.hpp"
#include "../../cmagi/gpsmoothing.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  const vec yobsInput = armaGetPrVec(prhs[0]);
  const mat distInput = armaGetPr(prhs[1]);

  char *str1;
  str1 = mxArrayToString(prhs[2]);
  string kernelInput = str1;

  const double sigmaExogenScalar = mxGetScalar(prhs[3]);
  bool useFrequencyBasedPrior = mxGetLogicals(prhs[4])[0];
  
  vec retv = gpsmooth(yobsInput, distInput, kernelInput, sigmaExogenScalar, useFrequencyBasedPrior);
  
  
  mwSize retdim[1];
  retdim[0] = retv.n_elem;
  
  plhs[0] = mxCreateNumericArray(1, retdim, mxDOUBLE_CLASS, mxREAL);
  
  armaSetPr(plhs[0], retv);
    

}