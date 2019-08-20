#include "mex.h"
#include "armaMex.hpp"
#include "../gpds_cpp/tgtdistr.h"
//#include "access.h"

using namespace arma;

// void matlab2arma(mat& A, const mxArray *mxdata){
//   access::rw(A.mem)=mxGetPr(mxdata);
//   access::rw(A.n_rows)=mxGetM(mxdata);
//   access::rw(A.n_cols)=mxGetN(mxdata);
//   access::rw(A.n_elem)=A.n_rows*A.n_cols;
// };
// 
// void matlab2armavec(vec& A, const mxArray *mxdata){
//   access::rw(A.mem)=mxGetPr(mxdata);
//   access::rw(A.n_cols)=mxGetN(mxdata);
//   access::rw(A.n_elem)=mxGetN(mxdata);
// }
// 
// void freeVar(mat& A, const double *ptr){
//     access::rw(A.mem)=ptr;
//     access::rw(A.n_rows)=1;
//     access::rw(A.n_cols)=1;
//     access::rw(A.n_elem)=1;
// };

// void armaSetPr(mxArray *matlabMatrix, const Mat<double>& armaMatrix) {
//         double *dst_pointer = mxGetPr(matlabMatrix);
//   const double *src_pointer = armaMatrix.memptr();
// 
//   std::memcpy(dst_pointer, src_pointer, sizeof(double)*armaMatrix.n_elem);
// }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 
//  vec phisig(1);
//  const double* phisigmem=access::rw(phisig.mem);
//  matlab2armavec(phisig, prhs[0]);
 
 vec phisig = armaGetPrVec(prhs[0]);
 //vec phisig(phisigmat);
 //phisig.print();
 
//  mat yobs(1,1);
//  const double* yobsmem=access::rw(yobs.mem);
//  matlab2arma(yobs, prhs[1]);
 mat yobs = armaGetPr(prhs[1]);
 //yobs.print();
 
 
 //mat dist(1,1);
 //const double* distmem=access::rw(dist.mem);
 //matlab2arma(dist, prhs[2]);
 mat dist = armaGetPr(prhs[2]);
 //dist.print();
 
 char *str1;
 str1 = mxArrayToString(prhs[3]);
 string kernel = str1;
 //mexPrintf("%s", str1);
 
     if(str1 == NULL )
        mexErrMsgTxt("Could not convert input to string.");   
 
 //cout << kernel;//print();
 
 //gpcov res = maternCov(phisig, dist,1);
//    int n = yobs.n_rows;
//   unsigned int obsDimension = yobs.n_cols;
//   int phiDimension = (phisig.size() - 1) / obsDimension;
//   double sigma = phisig(phisig.size() - 1);
//   const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()), 
//                               phiDimension, obsDimension, false, false);
//   
//   cout << n << " " << obsDimension << " "<< phiDimension << " "<< sigma << " " << phiAllDim;
//  // likelihood value part
//   std::function<gpcov(vec, mat, int)> kernelCov;
//   if(kernel == "matern"){
//     kernelCov = maternCov;
//   }else if(kernel == "rbf"){//  mat yobs(1,1);
//  const double* yobsmem=access::rw(yobs.mem);
//  matlab2arma(yobs, prhs[1]);

//     kernelCov = rbfCov;
//   }else if(kernel == "compact1"){
//     kernelCov = compact1Cov;
//   }else if(kernel == "generalMatern"){
//     kernelCov = generalMaternCov;
//   }else{
//     throw std::runtime_error("kernel is not specified correctly");
//   }
//   
//   lp ret;  
//   ret.gradient = zeros( phisig.size());
//   ret.value = 0;
//   
//   // V 
//   vec eigval(n);
//   mat eigvec(n, n);
//   
//   for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
//     gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
//     covThisDim.C.diag() += pow(sigma, 2);
//     
//     covThisDim.C.print();
//      eig_sym( eigval, eigvec, covThisDim.C );
//     vec eta = eigvec.t() * yobs.col(pDimEach);
//     ret.value += -n/2.0*log(2.0*datum::pi) - sum(log(eigval))/2.0 - 0.5*sum(square(eta) / eigval);
//     
//     vec alpha = eigvec * (eta / eigval);
//     mat facVtemp = alpha * alpha.t() - (eigvec.each_row() % (1.0 / eigval).t()) * eigvec.t();
//     // mat facVtemp = alpha * alpha.t() - inv( covThisDim.C);
//     double dVdsig = sigma * sum(facVtemp.diag());
//     vec dVdphi(covThisDim.dCdphiCube.n_slices);
//     for(unsigned int i=0; i < dVdphi.size(); i++){
//       ret.gradient(pDimEach*phiDimension + i) = accu(facVtemp % covThisDim.dCdphiCube.slice(i))/2.0;
//     }
//     ret.gradient(ret.gradient.size()-1) += dVdsig;
//   }  
//   
  //cout << "ok up to here\n";
 lp res = phisigllik(phisig, yobs, dist, kernel);
// 
//  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT16_CLASS, mxREAL);
//      int* data = (int*) mxGetData(plhs[0]); 
//     // This asigns data the same memory address as plhs[0]. the return value
// 
//     data[0]=123; //or data[0]=anyFunctionThatReturnsInt();
//     return;

 //plhs[0] = mxCreateDoubleMatrix(dist.n_rows, dist.n_cols, mxREAL); 
 //armaSetPr(plhs[0], dist);
 plhs[0] = mxCreateDoubleScalar(res.value);
 
 plhs[1] = mxCreateDoubleMatrix(res.gradient.n_rows, res.gradient.n_cols, mxREAL); 
 armaSetPr(plhs[1], res.gradient);
 
 
 //freeVar(phisig,phisigmem);
 //freeVar(dist,distmem);
 //freeVar(yobs,yobsmem);
 
 
}