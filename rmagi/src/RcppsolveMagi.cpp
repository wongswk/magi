#include "RcppArmadillo.h"
#include "classDefinition.h"
#include "Sampler.h"
#include "MagiSolver.h"
#include "dynamicalSystemModels.h"
#include "RcppTranslation.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List solveMagiRcpp(
        const arma::mat & yFull,
        const List & odeModel,
        const arma::vec & tvecFull,
        const arma::vec & sigmaExogenous,
        const arma::mat & phiExogenous,
        const arma::mat & xInitExogenous,
        const arma::vec & thetaInitExogenous,
        const arma::mat & muExogenous,
        const arma::mat & dotmuExogenous,
        const double priorTemperatureLevel,
        const double priorTemperatureDeriv,
        const double priorTemperatureObs,
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
    const Rcpp::String modelName = odeModel.containsElementNamed("name") ? as<Rcpp::String>(odeModel["name"]) : "";
    if(modelName == "FN"){
        modelC = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, arma::zeros(3), arma::ones(3)*INFINITY);
    }else if(modelName == "Hes1"){
        modelC = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, arma::zeros(7), arma::ones(7)*INFINITY);
    }else if(modelName == "Hes1-log"){
        modelC = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, arma::zeros(7), arma::ones(7)*INFINITY);
    }else if(modelName == "HIV"){
        modelC = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-10, 0,0,0,0,0, -10,-10,-10}, arma::ones(9)*10);
    }else if(modelName == "PTrans"){
        modelC = OdeSystem(ptransmodelODE, ptransmodelDx, ptransmodelDtheta, arma::zeros(6), arma::ones(6)*4);
    }else if(modelName == "Hes1-log-fixg"){
        modelC = OdeSystem(hes1logmodelODEfixg, hes1logmodelDxfixg, hes1logmodelDthetafixg, arma::zeros(6), arma::ones(6)*INFINITY);
    }else if(modelName == "Hes1-log-fixf"){
        modelC = OdeSystem(hes1logmodelODEfixf, hes1logmodelDxfixf, hes1logmodelDthetafixf, arma::zeros(6), arma::ones(6)*INFINITY);
    }else if(modelName == "Michaelis-Menten"){
        modelC = OdeSystem(MichaelisMentenModelODE, MichaelisMentenModelDx, MichaelisMentenModelDtheta, arma::zeros(3), arma::ones(3)*INFINITY);
    }else if(modelName == "Michaelis-Menten-Va"){
        modelC = OdeSystem(MichaelisMentenModelVaODE, MichaelisMentenModelVaDx, MichaelisMentenModelVaDtheta, arma::zeros(3), arma::ones(3)*INFINITY);
    }else if(modelName == "Michaelis-Menten-Vb4p"){
        modelC = OdeSystem(MichaelisMentenModelVb4pODE, MichaelisMentenModelVb4pDx, MichaelisMentenModelVb4pDtheta, arma::zeros(4), arma::ones(4)*INFINITY);
    }else if(modelName == "Michaelis-Menten-Vb2p"){
        modelC = OdeSystem(MichaelisMentenModelVb2pODE, MichaelisMentenModelVb2pDx, MichaelisMentenModelVb2pDtheta, arma::zeros(2), arma::ones(2)*INFINITY);
    }else if(modelName == "Michaelis-Menten-Log"){
        modelC = OdeSystem(MichaelisMentenLogModelODE, MichaelisMentenLogModelDx, MichaelisMentenLogModelDtheta, arma::zeros(3), arma::ones(3)*INFINITY);
    }else if(modelName == "lac-operon"){
        modelC = OdeSystem(lacOperonODE, lacOperonDx, lacOperonDtheta, arma::zeros(17), arma::ones(17)*INFINITY);
    }else if(modelName == "repressilator-gene-regulation"){
        modelC = OdeSystem(repressilatorGeneRegulationODE, repressilatorGeneRegulationDx, repressilatorGeneRegulationDtheta, arma::zeros(4), arma::ones(4)*INFINITY);
    }else if(modelName == "repressilator-gene-regulation-log"){
        modelC = OdeSystem(repressilatorGeneRegulationLogODE, repressilatorGeneRegulationLogDx, repressilatorGeneRegulationLogDtheta, arma::zeros(4), arma::ones(4)*INFINITY);
    }else{
        Rcpp::Rcout << "use R odeModel implementation\n";
        const Rcpp::Function & fOdeR = as<const Function>(odeModel["fOde"]);
        const Rcpp::Function & fOdeDxR = as<const Function>(odeModel["fOdeDx"]);
        const Rcpp::Function & fOdeDthetaR = as<const Function>(odeModel["fOdeDtheta"]);


        const Rcpp::NumericVector & thetaLowerBoundR = as<const NumericVector>(odeModel["thetaLowerBound"]);
        const Rcpp::NumericVector & thetaUpperBoundR = as<const NumericVector>(odeModel["thetaUpperBound"]);

        modelC.thetaUpperBound = arma::vec(const_cast<double*>( &(thetaUpperBoundR[0])), thetaUpperBoundR.size(), false, false);
        modelC.thetaLowerBound = arma::vec(const_cast<double*>( &(thetaLowerBoundR[0])), thetaLowerBoundR.size(), false, false);
        modelC.thetaSize = modelC.thetaLowerBound.size();

        modelC.fOde = [fOdeR](const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) -> arma::mat {
            return r2armamat(fOdeR(theta, x, tvec));
        };

        modelC.fOdeDx = [fOdeDxR](const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) -> arma::cube {
            return r2armacube(fOdeDxR(theta, x, tvec));
        };

        modelC.fOdeDtheta = [fOdeDthetaR](const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) -> arma::cube {
            return r2armacube(fOdeDthetaR(theta, x, tvec));
        };
    }

    MagiSolver solver(yFull,
                      modelC,
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
        Rcpp::Rcout << "phi = \n" << solver.phiAllDimensions << "\n";
    }
    solver.initXmudotmu();

    solver.initTheta();
    if(verbose){
        Rcpp::Rcout << "thetaInit = \n" << solver.thetaInit << "\n";
    }

    //return yFull;
    solver.initMissingComponent();
    solver.sampleInEpochs();

    Rcpp::List ret = Rcpp::List::create(Rcpp::Named("llikxthetasigmaSamples")=solver.llikxthetasigmaSamples,
                                        Rcpp::Named("phi")=solver.phiAllDimensions);

    return ret;

}
