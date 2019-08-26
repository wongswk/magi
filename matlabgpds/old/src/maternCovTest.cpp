#include "mex.h"
#include "tgtdistr.h"

using namespace arma;

// Helper stuff for Matlab-ARMA type conversions
void matlab2arma(mat& A, const mxArray *mxdata){
    access::rw(A.mem)=mxGetPr(mxdata);
    access::rw(A.n_rows)=mxGetM(mxdata);
    access::rw(A.n_cols)=mxGetN(mxdata);
    access::rw(A.n_elem)=A.n_rows*A.n_cols;
};

void matlab2armavec(vec& A, const mxArray *mxdata){
    access::rw(A.mem)=mxGetPr(mxdata);
    access::rw(A.n_cols)=mxGetN(mxdata);
    access::rw(A.n_elem)=mxGetN(mxdata);
}

void freeVar(mat& A, const double *ptr){
    access::rw(A.mem)=ptr;
    access::rw(A.n_rows)=1;
    access::rw(A.n_cols)=1;
    access::rw(A.n_elem)=1;
};

void armaSetPr(mxArray *matlabMatrix, const Mat<double>& armaMatrix) {
    double *dst_pointer = mxGetPr(matlabMatrix);
    const double *src_pointer = armaMatrix.memptr();

    std::memcpy(dst_pointer, src_pointer, sizeof(double)*armaMatrix.n_elem);
}
// end helper stuff

// main function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    vec phi(1);
    const double* phimem=access::rw(phi.mem);
    matlab2armavec(phi, prhs[0]);
    //phi.print();

    mat dist(1,1);
    const double* distmem=access::rw(dist.mem);
    matlab2arma(dist, prhs[1]);
    //dist.print();

    gpcov out = maternCov(phi, dist,0);
    //out.C.print();

    plhs[0] = mxCreateDoubleMatrix(out.C.n_rows, out.C.n_cols, mxREAL);
    armaSetPr(plhs[0], out.C);
    freeVar(phi,phimem);
    freeVar(dist,distmem);

    return;
}
