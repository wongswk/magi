#include "RcppArmadillo.h"
#include "classDefinition.h"
#include "Sampler.h"
#include "GpdsSolver.h"

using namespace Rcpp;

arma::mat r2armamat(const SEXP & x){
    const Rcpp::NumericMatrix & xtmp = as<const NumericMatrix>(x);
    return arma::mat(const_cast<double*>( xtmp.begin()), xtmp.nrow(), xtmp.ncol(), false, false);
}

arma::cube r2armacube(const SEXP & x){
    const Rcpp::NumericVector & xtmp = as<const NumericVector>(x);
    IntegerVector dim = xtmp.attr("dim");
    return arma::cube(const_cast<double*>( xtmp.begin()), dim[0], dim[1], dim[2], false, false);
}

// [[Rcpp::export]]
arma::mat solveGpdsRcpp(
        const arma::mat & yFull,
        const List & odeModel,
        const arma::vec & tvecFull,
        const arma::vec & sigmaExogenous,
        const arma::mat & phiExogenous,
        const double priorTemperatureLevel,
        const double priorTemperatureDeriv,
        const std::string kernel,
        const int nstepsHmc,
        const double burninRatioHmc,
        const unsigned int niterHmc,
        const double stepSizeFactorHmc,
        const int nEpoch,
        const int bandSize,
        const bool useFrequencyBasedPrior,
        const bool useBand,
        const bool useMean,
        const bool useScalerSigma,
        const bool useFixedSigma,
        const bool verbose) {

    OdeSystem modelC;
    //model = OdeSystem(fOde, fOdeDx, fOdeDtheta, lb, ub);
    //const Rcpp::Function fOdeR = odeModel["fOde"];
    //const Rcpp::Function fOdeDxR = odeModel["fOdeDx"];
    //const Rcpp::Function fOdeDthetaR = odeModel["fOdeDtheta"];
    const Rcpp::Function & fOdeR = as<const Function>(odeModel["fOde"]);
    const Rcpp::Function & fOdeDxR = as<const Function>(odeModel["fOdeDx"]);
    const Rcpp::Function & fOdeDthetaR = as<const Function>(odeModel["fOdeDtheta"]);


    const Rcpp::NumericVector & thetaLowerBoundR = as<const NumericVector>(odeModel["thetaLowerBound"]);
    const Rcpp::NumericVector & thetaUpperBoundR = as<const NumericVector>(odeModel["thetaUpperBound"]);

    modelC.thetaUpperBound = arma::vec(const_cast<double*>( &(thetaUpperBoundR[0])), thetaUpperBoundR.size(), false, false);
    modelC.thetaLowerBound = arma::vec(const_cast<double*>( &(thetaLowerBoundR[0])), thetaLowerBoundR.size(), false, false);
    modelC.thetaSize = modelC.thetaLowerBound.size();

    modelC.fOde = [& fOdeR](const arma::vec & theta, const arma::mat & x) -> arma::mat {
        return r2armamat(fOdeR(theta, x));
    };

    modelC.fOdeDx = [& fOdeDxR](const arma::vec & theta, const arma::mat & x) -> arma::cube {
        return r2armacube(fOdeDxR(theta, x));
    };

    modelC.fOdeDtheta = [& fOdeDthetaR](const arma::vec & theta, const arma::mat & x) -> arma::cube {
        return r2armacube(fOdeDthetaR(theta, x));
    };

    GpdsSolver solver(yFull,
                      modelC,
                      tvecFull,
                      sigmaExogenous,
                      phiExogenous,
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

    //return yFull;
    solver.initMissingComponent();
    solver.doHMC();
    return arma::join_vert(solver.llikSamples.t(), solver.xthetasigmaSamples);

}


