#include "classDefinition.h"

class Sampler {
    const arma::mat & yobs;
    const std::vector<gpcov> & covAllDimensions;
    const int nsteps;
    const bool traj;
    const std::string loglikflag;
    const arma::vec priorTemperature;
    const OdeSystem & model;
    const unsigned int sigmaSize;
    const double burninRatio;
    const double niter;
    bool useBand;
    bool useMean;
    std::function<lp(arma::vec)> tgt;
    arma::vec lb, ub;
    arma::vec lliklist;
    arma::mat xth;
public:
    hmcstate sampleSingle(const arma::vec &xthetasigmaInit, const arma::vec & step);
    void sampleChian(const arma::vec &xthetasigmaInit, const arma::vec &stepLowInit, bool verbose);
    Sampler(const arma::mat & yobsInput,
            const std::vector<gpcov> & covAllDimensionsInput,
            const int nstepsInput,
            const std::string loglikflagInput,
            const arma::vec priorTemperatureInput,
            const OdeSystem & modelInput,
            const unsigned int niterInput,
            const double burninRatioInput);
};
