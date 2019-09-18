# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>

"""
After experiments have been created using the script createExperiments.py and
the hyperparameters were trained using getHyperparams.py,
this script can be run to infer the ODE parameters using FGPGM.
"""
plotting=True     #
useMagnitudes = False 

import numpy as np

from FGPGM.Experiments.FitzHughNagumo import FHN as exp
from FGPGM.Kernels.Matern52 import Matern52 as kernel

from FGPGM.FGPGM import FGPGM

"""determine standardization"""
standardize = np.loadtxt("standardize.csv")
if standardize == 0:
    standardize=False
    print("Do not use standardization of the observations")
elif standardize == 1:
    standardize=True
    print("Use standardization of the observations")
else:
    raise ValueError("Illegal value encountered in standardize.csv." +
                     "Maybe the file is corrupted?")

"""experiment and observation loading"""
time = np.loadtxt("time.csv")

# create experiment
experiment = exp()
# load noisy states at observation times
y = np.loadtxt("observations.csv")
# load hyperparameter
gammaValue = np.loadtxt("gammaValues.csv") # state derivative in the GP model. 1e-4 ood value
gammas = gammaValue*np.ones(y.shape[1])

"""create kernels using hyperparameters calculated by getHyperparams.py"""
kernels = []
for state in np.arange(y.shape[1]):
    currentKernel = kernel(theta=np.abs(np.random.randn(3)),
                           sigma=np.abs(np.random.randn(1)))
    currentKernel.theta = np.loadtxt("hyperparams/theta{}.csv".format(state))
    currentKernel.sigma = np.squeeze(np.loadtxt(
        "hyperparams/sigma{}.csv".format(state))) 
    kernels.append(currentKernel)        

"""find initial values for theta"""
trueTheta = np.loadtxt("trueODEParams.csv")
theta0 = np.abs(np.random.randn(trueTheta.size))    

"""create FGPGM object"""
FM = FGPGM(kernels=kernels,
           time=time,
           y=y,
           experiment=experiment,
           nODEParams=trueTheta.size,
           gamma=gammas,
           normalize=True,
           standardize=standardize)

"""perform actual parameter estimation"""
stateStds = np.ones(y.size)*0.075
paramStds = np.ones(theta0.size)*0.12
propStds = np.concatenate([stateStds, paramStds])

newStates, newParams = FM.getFGPGMResults(GPPosteriorInit=True,
                                          blockNegStates=False,
                                          debug=True,
                                          theta0=theta0,
                                          thetaMagnitudes=np.zeros_like(theta0),
                                          nSamples=300000,
                                          nBurnin=1000,
                                          propStds=propStds)
np.savetxt("optimalParamsFGPGM.csv", newParams)
np.savetxt("optimalStatesFGPGM.csv", newStates)
print(newParams)
