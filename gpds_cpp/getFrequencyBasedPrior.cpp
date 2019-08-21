#include <armadillo>

arma::vec getFrequencyBasedPrior(const arma::vec & x){
    arma::cx_mat z = arma::fft(x);
    arma::vec zmod(z.n_rows);
    for(unsigned int i = 0; i < z.n_rows; i++){
        zmod(i) = std::norm(z(i, 0));
    }
    const arma::vec & zmodEffective = zmod(arma::span(1, (zmod.size() - 1) / 2));
    const arma::vec zmodEffectiveSorted = arma::sort(zmodEffective);
    double upperQuarter = zmodEffectiveSorted(static_cast<unsigned int>(std::ceil(zmodEffectiveSorted.size() * 0.75)));
    double lowerQuarter = zmodEffectiveSorted(static_cast<unsigned int>(std::floor(zmodEffectiveSorted.size() * 0.25)));
    double iqr = upperQuarter - lowerQuarter;
    arma::vec outliers = zmodEffective(zmodEffective > upperQuarter + 1.5 * iqr);
    long long int freq = 1;
    if(!outliers.empty()){
        for(unsigned long long int i = zmodEffective.size(); i > 0; i--){
            if(arma::min(arma::abs(zmodEffective(i - 1) - outliers)) < 1e-6){
                freq = i;
                break;
            }
        }

    }else{
        freq = zmodEffective.index_max();
    }

    double meanFactor = 0.5 / freq;
    double sdFactor = (1 - meanFactor) / 3;

    return arma::vec({meanFactor, sdFactor});
}
