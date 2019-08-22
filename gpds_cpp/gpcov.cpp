#include <armadillo>

#include "tgtdistr.h"
#include "classDefinition.h"

arma::mat mat2band(const arma::mat & matInput, const int bandsize){
    int ndim = matInput.n_rows;
    arma::mat matOutput(2 * bandsize + 1, ndim);
    for (int j = 1; j <= matInput.n_cols; j++){
        int k = bandsize + 1 - j;
        for (int i = std::max(1, j - bandsize); i <= std::min(ndim, j + bandsize); i++){
            matOutput(k+i-1, j-1) = matInput(i-1, j-1);
        }
    }
    return matOutput;
}

void gpcov::addBandCov(const int bandsizeInput){
    bandsize = bandsizeInput;


}
