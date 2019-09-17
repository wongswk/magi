# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the hyperparameters needed by doFGPGM.py. Should only
be run after the experiment has been created using createExperiments.py.
"""
import numpy as np
from FGPGM.Kernels.Sigmoid import Sigmoid
import os

from matplotlib import pyplot as plt
plt.switch_backend('agg')

# decide if the data should be standardized in a preprocessing step
standardize = True
np.savetxt("standardize.csv", np.asarray(1).reshape([1, 1]))

"""kernel optimization"""
# amount of basinhopping iterations in kernel hyperparameter optimization
kernelIter = 100  # should be adjusted based on roughness of likelihood
gammaValue = 1e-4
np.savetxt("gammaValues.csv", np.asarray(gammaValue).reshape([-1, 1]))
# data from createExperiment.py
y = np.loadtxt("observations.csv")
time = np.loadtxt("time.csv")

"""create Kernel and Kernel matrices"""
for state in np.arange(y.shape[1]):
    currentKernel = Sigmoid(theta=np.abs(np.random.randn(3)),
                        sigma=np.abs(np.random.randn(1)))
    currentKernel.learnHyperparams(
        currentKernel.theta,
        currentKernel.sigma,
        y[:, state],
        time,
        normalize=True,
        standardize=standardize,
        newNugget=1e-4,
        anneal=False,
        basinIter=kernelIter
        )
# save results for later use by doFGPGM.py
    if not os.path.exists("./hyperparams"):
        os.makedirs("./hyperparams")
    np.savetxt("hyperparams/theta{}.csv".format(state),
               currentKernel.theta)
    np.savetxt("hyperparams/sigma{}.csv".format(state),
               np.asarray(currentKernel.sigma).reshape([1, 1]))
