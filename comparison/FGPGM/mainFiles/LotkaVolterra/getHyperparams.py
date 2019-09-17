# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the hyperparameters needed by doFGPGM.py. Should only
be run after the experiment has been created using createExperiments.py.
"""
import numpy as np
from FGPGM.Kernels.RBF import RBF
import os

from matplotlib import pyplot as plt
plt.switch_backend('agg')

# decide if the data should be standardized in a preprocessing step
standardize = True
np.savetxt("standardize.csv", np.asarray(1).reshape([1, 1]))

"""kernel optimization"""
# amount of basinhopping iterations in kernel hyperparameter optimization
kernelIter = 100  # should be adjusted based on roughness of likelihood
gammaValue = 3e-1
np.savetxt("gammaValues.csv", np.asarray(gammaValue).reshape([-1, 1]))
# data from createExperiment.py
y = np.loadtxt("observations.csv")
time = np.loadtxt("time.csv")

"""create Kernel and Kernel matrices"""
kernel = RBF(theta=np.abs(np.random.randn(2)),
             sigma=np.abs(np.random.randn(1)))
kernel.learnHyperparams(
    kernel.theta,
    kernel.sigma,
    y,
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
np.savetxt("hyperparams/theta.csv",
           kernel.theta)
np.savetxt("hyperparams/sigma.csv",
           np.asarray(kernel.sigma).reshape([1, 1]))
