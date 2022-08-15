#include "tgtdistr.h"
#include "band.h"
#include "dynamicalSystemModels.h"
#include <boost/math/special_functions/bessel.hpp>

using namespace arma;

//' leave one out cross validation for Gaussian Process fitting with various kernel
//'
//' loss function is predictive log likelihood
//'
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigloocvllik( const vec & phisig,
                    const mat & yobs,
                    const mat & dist,
                    string kernel){
    int n = yobs.n_rows;
    unsigned int obsDimension = yobs.n_cols;
    int phiDimension = (phisig.size() - 1) / obsDimension;
    const double sigma = phisig(phisig.size() - 1);
    const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()),
                                phiDimension, obsDimension, true, false);

    std::function<gpcov(vec, mat, int)> kernelCov;
    if(kernel == "matern"){
        kernelCov = maternCov;
    }else if(kernel == "rbf"){
        kernelCov = rbfCov;
    }else if(kernel == "compact1"){
        kernelCov = compact1Cov;
    }else if(kernel == "periodicMatern"){
        kernelCov = periodicMaternCov;
    }else if(kernel == "generalMatern"){
        kernelCov = generalMaternCov;
    }else{
        throw std::runtime_error("kernel is not specified correctly");
    }

    lp ret;
    ret.gradient = zeros( phisig.size());
    ret.value = 0;

    for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
        gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
        covThisDim.C.diag() += pow(sigma, 2);

        mat Cinv = arma::inv_sympd(covThisDim.C);

        vec alpha = Cinv * yobs.col(pDimEach);

        vec muLooCv = yobs.col(pDimEach) - alpha / Cinv.diag();
        vec sigmasqLooCv = 1 / Cinv.diag();

        ret.value += -n/2.0*log(2.0*datum::pi) - 0.5*sum(log(sigmasqLooCv)) -
                     0.5*sum(square(yobs.col(pDimEach) - muLooCv) / sigmasqLooCv);

        cube Ztemp = Cinv * covThisDim.dCdphiCube.each_slice();
        Ztemp = join_slices(Ztemp, 2*sigma*Cinv);
        for(unsigned int j=0; j < Ztemp.n_slices; j++){
            vec ZjCinvDiag = sum(Ztemp.slice(j).t() % Cinv).t();
            double dLdphisig = sum((alpha % (Ztemp.slice(j) * alpha) - 0.5*(1+square(alpha)/Cinv.diag())%ZjCinvDiag)/Cinv.diag());
            if(j < covThisDim.dCdphiCube.n_slices){
                ret.gradient(pDimEach*phiDimension + j) = dLdphisig;
            }else{
                ret.gradient(ret.gradient.size()-1) += dLdphisig;
            }
        }
    }
    return ret;
}

//' leave one out cross validation for Gaussian Process fitting with various kernel
//'
//' loss function is predictive log likelihood
//'
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigloocvmse( const vec & phisig,
                   const mat & yobs,
                   const mat & dist,
                   string kernel){
    unsigned int obsDimension = yobs.n_cols;
    int phiDimension = (phisig.size() - 1) / obsDimension;
    const double sigma = phisig(phisig.size() - 1);
    const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()),
                                phiDimension, obsDimension, true, false);

    std::function<gpcov(vec, mat, int)> kernelCov;
    if(kernel == "matern"){
        kernelCov = maternCov;
    }else if(kernel == "rbf"){
        kernelCov = rbfCov;
    }else if(kernel == "compact1"){
        kernelCov = compact1Cov;
    }else if(kernel == "periodicMatern"){
        kernelCov = periodicMaternCov;
    }else if(kernel == "generalMatern"){
        kernelCov = generalMaternCov;
    }else{
        throw std::runtime_error("kernel is not specified correctly");
    }

    lp ret;
    ret.gradient = zeros( phisig.size());
    ret.value = 0;

    for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
        gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
        covThisDim.C.diag() += pow(sigma, 2);

        mat Cinv = arma::inv_sympd(covThisDim.C);

        vec alpha = Cinv * yobs.col(pDimEach);

        vec muLooCv = yobs.col(pDimEach) - alpha / Cinv.diag();
        vec sigmasqLooCv = 1 / Cinv.diag();

        ret.value += -0.5*sum(square(yobs.col(pDimEach) - muLooCv));

        cube Ztemp = Cinv * covThisDim.dCdphiCube.each_slice();
        Ztemp = join_slices(Ztemp, 2*sigma*Cinv);
        for(unsigned int j=0; j < Ztemp.n_slices; j++){
            vec ZjCinvDiag = sum(Ztemp.slice(j).t() % Cinv).t();
            vec dmudphisig = (Ztemp.slice(j) * alpha) / Cinv.diag() - alpha % ZjCinvDiag / square(Cinv.diag());
            if(j < covThisDim.dCdphiCube.n_slices){
                ret.gradient(pDimEach*phiDimension + j) = sum(dmudphisig % (yobs.col(pDimEach) - muLooCv));
            }else{
                ret.gradient(ret.gradient.size()-1) += sum(dmudphisig % (yobs.col(pDimEach) - muLooCv));
            }
        }
    }
    return ret;
}
