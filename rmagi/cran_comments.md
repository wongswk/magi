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
