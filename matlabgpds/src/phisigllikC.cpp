#include "mex.h"
#include "armaMex.hpp"
#include "../gpds_cpp/tgtdistr.h"

using namespace arma;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 
 
 vec phisig = armaGetPrVec(prhs[0]);
 mat yobs = armaGetPr(prhs[1]);
 mat dist = armaGetPr(prhs[2]);
 
 char *str1;
 str1 = mxArrayToString(prhs[3]);
 string kernel = str1;
 
     if(str1 == NULL )
        mexErrMsgTxt("Could not convert input to string.");   
 
 lp res = phisigllik(phisig, yobs, dist, kernel);

 plhs[0] = mxCreateDoubleScalar(res.value);
 
 plhs[1] = mxCreateDoubleMatrix(res.gradient.n_rows, res.gradient.n_cols, mxREAL); 
 armaSetPr(plhs[1], res.gradient);
 
 
 
}