# MAGI v.1.2.3

## Resubmission
This is a resubmission to update citations and references. The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.

## Changelog:
* Add CITATION file and update references for a new JSS publication
* Minor bugfix for theta argument in MagiSolver

## Test environments:
* Ubuntu 20.04, R 4.4.0
* Windows 10, R 4.4.0
* Windows 10, R-devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was a NOTE regarding the DOI that will be registered after publication on CRAN.
```
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.18637/jss.v109.i04
    From: DESCRIPTION
    Status: 404
    Message: Not Found

Found the following (possibly) invalid DOIs:
  DOI: 10.18637/jss.v109.i04
    From: DESCRIPTION
          inst/CITATION
          man/MagiSolver.Rd
          man/magi.Rd
    Status: 404
    Message: Not Found
```

There was a NOTE on R 4.4.0 on Ubuntu 20.04 regarding package size:
```
* checking installed package size ... NOTE
  installed size is 54.4Mb
  sub-directories of 1Mb or more:
    libs  53.0Mb
```


# MAGI v.1.2.2

## Resubmission
This is a resubmission. In this version we have made the following changes.

## Changelog:
* MagiSolver supports more covariance kernels (rbf, compact1, matern, periodicMatern) that can be selected via the kerneltype option
* Added functionality to plot.magioutput to generate MCMC traceplots
* Allow selection of posterior mean/median/mode as point estimates in plot/summary functions
* Added helper functions to compute GP mean and covariance given observations for easier visualization of initial GP fit
* Removed specification of C++11 in Makevars
* Minor documentation updates/clarifications

## Test environments:
* Ubuntu 20.04, R 4.2.2
* Windows 10, R 4.3.0
* Windows 10, R-devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE on R 4.2.2 on Ubuntu 20.04 regarding package size:
```
* checking installed package size ... NOTE
  installed size is 40.3Mb
  sub-directories of 1Mb or more:
    libs  38.9Mb
```


# MAGI v1.2.1

## Resubmission
This is a resubmission. In this version we have made the following changes.

## Changelog:
* Options added to MagiSolver: vector input for HMC step sizes, non-negative systems, skip missing component initialization

## Test environments:
* Ubuntu 20.04, R 4.2.2
* Windows 10, R 4.2.2
* Windows 10, R-devel

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.

# MAGI v1.2.0

## Resubmission
This is a resubmission. In this version we have made the following changes.

## Changelog:
* MagiSolver now returns a classed object (magioutput) with associated print, summary, plot methods
* Added .Rd documentation for three ODE system examples
* Added an example dataset based on Fitzhugh-Nagumo equations
* Cleaned up unused functions and code from the package
* Added verbose option to reduce clutter in console when running MagiSolver
* Replaced solve by forwardsolve where possible for efficiency
* Reduced runtime of examples

## Test environments:
* Ubuntu 20.04, R 4.2.0
* Windows 10, R 4.1.0
* Windows 10, R-devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE on R 4.2.0 on Ubuntu 20.04 regarding package size:
```
* checking installed package size ... NOTE
  installed size is 38.6Mb
  sub-directories of 1Mb or more:
    libs  37.2Mb
```

There was one NOTE on R 4.1.0 on Windows 10 regarding package size:
```
> checking installed package size ... NOTE
    installed size is  5.8Mb
    sub-directories of 1Mb or more:
      libs   4.5Mb
```

There were no NOTEs on R-devel on Windows 10.
