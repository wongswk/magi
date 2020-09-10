# dynamic-systems
[![Travis-CI Build Status](https://travis-ci.com/Shihao-Yang/dynamic-systems.svg?token=zsECgNMyrthwbokp6yPB&branch=master)](https://travis-ci.com/Shihao-Yang/dynamic-systems)
[![codecov](https://codecov.io/gh/Shihao-Yang/dynamic-systems/branch/master/graph/badge.svg?token=Sr7hFVaajH)](https://codecov.io/gh/Shihao-Yang/dynamic-systems)

This repository contains the accompanying software for the paper "Inference of dynamic systems via constrained Gaussian processes" by Shihao Yang, Samuel W.K. Wong, and S. C. Kou.

## Installation

User interfaces are available in R, MATLAB, and Python.

A comprehensive shell script `build.sh` is provided, which by default prepares all three interfaces.  Edit `build.sh` to specify the location of your R libraries, and remove the compilation blocks for any of R, MATLAB, Python that either will not be used or is unavailable on your system.  Then execute `build.sh` to install dependencies and compile the library.

The pre-compiled binary for c++, R, and python is available as a docker image on docker hub: https://hub.docker.com/repository/docker/shihaoyangphd/magi

## Usage

Inference is performed via the unified function `MagiSolver` which can be called from R, MATLAB, Python.  A description of its basic syntax is as follows, where D is the number of components in the dynamic system, and |I| is the number of discretization points for computation.

     MagiSolver(
       yFull,                 # |I|-by-D data matrix of observations Y, with unobserved entries set to NA
       model,                 # ODE model specification (see examples)
       tvecFull,              # length |I| vector of time points
       sigmaExogenous,        # (optional) length D vector of starting values of Gaussian noise SD sigma,
                                   recommended value: supply if known, else leave blank
       phiExogenous,          # (optional) 2-by-D matrix of GP hyperparameters phi
                                   recommended value: supply if known, else leave blank
       xInitExogenous,        # (optional) |I|-by-D matrix of starting values for X_I
                                   recommended value: linearly interpolate between observed points of Y
       thetaInitExogenous,    # (optional) starting value of theta
                                   recommended value: leave blank
       muExogenous,           # (optional) starting GP mean curve
                                   recommended value: leave blank
       dotmuExogenous,        # (optional) starting GP derivative of mean curve
                                   recommended value: leave blank
       priorTemperatureLevel, # tempering factor on GP prior
                                   recommended value: D|I| / (number of observed data values in Y)
       priorTemperatureDeriv, # tempering factor on GP derivative
                                   recommended value: D|I| / (number of observed data values in Y)
       priorTemperatureObs,   # tempering factor on observations
                                   recommended value: 1
       kernel,                # currently supported GP kernel is "generalMatern"
       nstepsHmc,             # number of leapfrog steps per HMC iteration
       burninRatioHmc,        # proportion of HMC iterations to treat as burn-in
       niterHmc,              # number of HMC iterations to run
       stepSizeFactorHmc,     # initial step size for HMC sampler
       nEpoch,                # currently supported value is 1
       bandSize,              # band size for band matrix approximation
       useFrequencyBasedPrior,  # recommended value: TRUE
       useBand,               # recommended value: TRUE
       useMean,               # recommended value: TRUE
       useScalerSigma,        # set to FALSE each component has its own noise level sigma, else TRUE
       useFixedSigma,         # set to TRUE if sigma known (must supply sigmaExogenous)
       verbose)               # set to TRUE to print out additional diagnostic information



## Examples

See the README in the corresponding subfolders: `rmagi` (for R), `pymagi` (for Python), `matlabmagi` (for Matlab).

There, we provide specific examples of how to set up and call `MagiSolver` in each software environment, and how to supply your own ODE systems and data to the method.

### Reference

For a full discussion of our method and examples, please see our paper "Inference of dynamic systems via constrained Gaussian processes" (https://arxiv.org/abs/xxxxx). 
