#include "RcppArmadillo.h"
#include "classDefinition.h"

namespace Rcpp
{
    // Date
    template <>
    gpcov as(SEXP x)
    {
        List cov_r(x);

        gpcov cov_v;

        // std::cout << "extract list as constant reference\n";
        const Rcpp::NumericMatrix & Cinv = as<const NumericMatrix>(cov_r["Cinv"]);
        const Rcpp::NumericMatrix & mphi = as<const NumericMatrix>(cov_r["mphi"]);
        const Rcpp::NumericMatrix & Kinv = as<const NumericMatrix>(cov_r["Kinv"]);
        const Rcpp::NumericMatrix & CinvBand = as<const NumericMatrix>(cov_r["CinvBand"]);
        const Rcpp::NumericMatrix & mphiBand = as<const NumericMatrix>(cov_r["mphiBand"]);
        const Rcpp::NumericMatrix & KinvBand = as<const NumericMatrix>(cov_r["KinvBand"]);
        const Rcpp::NumericVector & mu = as<const NumericVector>(cov_r["mu"]);
        const Rcpp::NumericVector & dotmu = as<const NumericVector>(cov_r["dotmu"]);

        // *(const_cast<double*>( &(mu[0]))) = -1; // this part is working -- R value changed

        // std::cout << "use R memory without copy\n";
        cov_v.Cinv = arma::mat(const_cast<double*>( Cinv.begin()), Cinv.nrow(), Cinv.ncol(), false, false);
        cov_v.mphi = arma::mat(const_cast<double*>( mphi.begin()), mphi.nrow(), mphi.ncol(), false, false);
        cov_v.Kinv = arma::mat(const_cast<double*>( Kinv.begin()), Kinv.nrow(), Kinv.ncol(), false, false);
        cov_v.CinvBand = arma::mat(const_cast<double*>( CinvBand.begin()), CinvBand.nrow(), CinvBand.ncol(), false, false);
        cov_v.mphiBand = arma::mat(const_cast<double*>( mphiBand.begin()), mphiBand.nrow(), mphiBand.ncol(), false, false);
        cov_v.KinvBand = arma::mat(const_cast<double*>( KinvBand.begin()), KinvBand.nrow(), KinvBand.ncol(), false, false);
        cov_v.mu = arma::vec(const_cast<double*>( &(mu[0])), mu.size(), false, false);
        cov_v.dotmu = arma::vec(const_cast<double*>( &(dotmu[0])), dotmu.size(), false, false);
        cov_v.bandsize = as<int>(cov_r["bandsize"]);

        // cov_v.mu(1) = 2; // this part is also working -- R value changed
        // cov_v.CinvBand(0) = 999;
        return cov_v;
    }

    template <>
    SEXP wrap(const gpcov& object)
    {
        return List::create(
            Named("C")=object.C,
            Named("dCdphiCube")=object.dCdphiCube,
            Named("Cprime")=object.Cprime,
            Named("Cdoubleprime")=object.Cdoubleprime,
            Named("dCprimedphiCube")=object.dCprimedphiCube,
            Named("dCdoubleprimedphiCube")=object.dCdoubleprimedphiCube,
            Named("Cinv")=object.Cinv,
            Named("mphi")=object.mphi,
            Named("Kinv")=object.Kinv,
            Named("Sigma")=object.Sigma,
            Named("dSigmadphiCube")=object.dSigmadphiCube
        );
    }
}