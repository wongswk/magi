# MAnifold-constrained Gaussian process Inference (MAGI)
[![Travis-CI Build Status](https://travis-ci.com/Shihao-Yang/dynamic-systems.svg?token=zsECgNMyrthwbokp6yPB&branch=master)](https://travis-ci.com/Shihao-Yang/dynamic-systems)
[![codecov](https://codecov.io/gh/Shihao-Yang/dynamic-systems/branch/master/graph/badge.svg?token=Sr7hFVaajH)](https://codecov.io/gh/Shihao-Yang/dynamic-systems)

This repository contains the accompanying software for the paper "[Inference of dynamic systems from noisy and sparse data via manifold-constrained Gaussian processes](https://doi.org/10.1073/pnas.2020397118)" by Shihao Yang, Samuel W.K. Wong, and S. C. Kou.

## Installation

User interfaces are available in R, MATLAB, and Python.

*For R users*: the current version of MAGI is now available as a package on CRAN, and may be installed via

`install.packages("magi")`

Please see https://cran.r-project.org/package=magi for the full documentation and vignette of this latest version, which contains a number of enhancements to facilitate ease of usage.

*For MATLAB users*: Run `build.sh` to first build the C++ library. Then navigate to `matlabmagi` and run `matlab_build.sh` to compile the required MEX files. See README in `matlabmagi` for further details.

*For Python users*: Run `build.sh` to first build the C++ library. Then navigate to `pymagi` and run `py_build.sh` to build the corresponding Python library.

Pre-compiled binaries for C++, R, and Python are also available as a Docker image on Docker Hub: https://hub.docker.com/repository/docker/shihaoyangphd/magi

## Usage

Inference is performed via the unified function `MagiSolver` which can be called from R, MATLAB, Python. The basic syntax is

```
MagiSolver(y, odeModel, control)
```
where `y` is the a data matrix, `odeModel` specifies the functions and parameters of the system, and `control` passes additional options.

For fully-documented examples and details, please refer to our MAGI software manuscript: https://arxiv.org/abs/2203.06066

## Examples

See the README in the corresponding subfolders: `rmagi` (for R), `pymagi` (for Python), `matlabmagi` (for Matlab).

There, we provide specific examples of how to set up and call `MagiSolver` in each software environment, and how to supply your own ODE systems and data to the method.

### Reference

For a full discussion of the method, please see our paper "Inference of dynamic systems from noisy and sparse data via manifold-constrained Gaussian processes", PNAS 118 (15), e2020397118 (https://doi.org/10.1073/pnas.2020397118).

For a full discussion of the software package and examples, please see our accompanying manuscript "MAGI: A Package for Inference of Dynamic Systems from Noisy and Sparse Data via Manifold-constrained Gaussian Processes", https://arxiv.org/abs/2203.06066
