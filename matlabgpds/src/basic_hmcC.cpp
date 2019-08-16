//
// Created by Shihao Yang on 2019-08-15.
//

#include "mex.h"
#include "matrix.h"

#include "hmc.h"
#include "classDefinition.h"
#include "../../rgpds/src/classDefinition.h"

using namespace arma;

// Helper stuff for Matlab-ARMA type conversions
void matlab2arma(mat& A, const mxArray *mxdata){
    access::rw(A.mem)=mxGetPr(mxdata);
    access::rw(A.n_rows)=mxGetM(mxdata);
    access::rw(A.n_cols)=mxGetN(mxdata);
    access::rw(A.n_elem)=A.n_rows*A.n_cols;
};

void matlab2armavec(vec& A, const mxArray *mxdata){
    A = vec(const_cast<double*>( mxGetPr(mxdata)), static_cast<int>(mxGetNumberOfElements(mxdata)), false, false);
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Argument Checking:
    // First Argument should be a Function Handle
    if( !mxIsClass( prhs[0] , "function_handle")) {
        mexErrMsgTxt("First input argument is not a function handle.");
    }
    // Second Argument is a Double Vector
    if (!mxIsClass(prhs[1], "double")||(mxGetM(prhs[1])>1)) {
        mexErrMsgTxt("Second input argument is not a double vector");
    }
    // 3rd Argument is a Double Vector
    if (!mxIsClass(prhs[2], "double")||(mxGetM(prhs[2])>1)) {
        mexErrMsgTxt("3rd input argument is not a double vector");
    }
    // 4th Argument is a Double Vector
    if (!mxIsClass(prhs[3], "double")||(mxGetM(prhs[3])>1)) {
        mexErrMsgTxt("4th input argument is not a double vector");
    }
    // 5th Argument is a Double Vector
    if (!mxIsClass(prhs[4], "double")||(mxGetM(prhs[4])>1)) {
        mexErrMsgTxt("5th input argument is not a double vector");
    }

    //processing on input arguments
    mxArray * lpr_matlab = const_cast<mxArray *>(prhs[0]);

    const std::function<lp (vec)> & lpr = [& lpr_matlab](const vec & x) -> lp {
        mxArray * x_matlab = mxCreateDoubleMatrix(x.size(), 1, mxREAL);
        armaSetPr(x_matlab, x);

        mxArray *rhs_lpr[2];
        rhs_lpr[0] = lpr_matlab;
        rhs_lpr[1] = x_matlab;

        mxArray *lhs_lpr[2];
        mexCallMATLAB(2,lhs_lpr,2,rhs_lpr,"feval");
        lp ret;
        ret.value = *mxGetPr(lhs_lpr[0]);
        matlab2armavec(ret.gradient, lhs_lpr[1]);

        mxDestroyArray(x_matlab);
        mxDestroyArray(lhs_lpr[0]);
        mxDestroyArray(lhs_lpr[1]);
        return ret;
    };

    vec initial, step, lb, ub;
    matlab2armavec(initial, prhs[1]);
    matlab2armavec(step, prhs[2]);
    matlab2armavec(lb, prhs[3]);
    matlab2armavec(ub, prhs[4]);
    int nsteps = static_cast<int>(*mxGetPr(prhs[5]));
    const bool traj = mxGetLogicals(prhs[6])[0];

    std::cout << "initial = " << initial << "\n";
    std::cout << "step = " << step << "\n";
    std::cout << "lb = " << lb << "\n";
    std::cout << "ub = " << ub << "\n";
    std::cout << "nsteps = " << nsteps << "\n";
    std::cout << "traj = " << traj << "\n";

    const hmcstate & hmc_out = basic_hmcC(lpr, initial, step, lb, ub, nsteps, traj);

    plhs[0] = mxCreateDoubleMatrix(hmc_out.final.size(), 1, mxREAL);
    armaSetPr(plhs[0], hmc_out.final);
    plhs[1] = mxCreateDoubleScalar(hmc_out.lprvalue);
}
