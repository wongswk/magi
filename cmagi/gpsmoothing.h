//
// Created by Shihao Yang on 2019-08-22.
//

#ifndef MAGI_MULTI_LANG_GPSMOOTHING_H
#define MAGI_MULTI_LANG_GPSMOOTHING_H
#define EIGEN_NO_DEBUG

arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   const double sigmaExogenScalar = -1,
                   bool useFrequencyBasedPrior = false);

#endif //MAGI_MULTI_LANG_GPSMOOTHING_H
