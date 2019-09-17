# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the hyperparameters used later FGPGM.
Must be run from the testRun{}.format{i} directories.
"""

import numpy as np
from FGPGM.Kernels.Matern52 import Matern52
import os

standardize = True

if standardize:
    np.savetxt("standardize.csv", np.asarray(1).reshape([1, 1]))
    print("standardizing the observations")
else:
    np.savetxt("standardize.csv", np.asarray(0).reshape([1, 1]))
    print("not standardizing the observations")

from matplotlib import pyplot as plt
plt.switch_backend('agg')

"""kernel decision"""
kernelIter = 100
gammaValue = np.loadtxt("gammaInput.csv")  #3e-1  # state derivative in the GP model.

np.savetxt("gammaValues.csv", np.asarray(gammaValue).reshape([-1, 1]))

y = np.loadtxt("observations.csv")
time = np.loadtxt("time.csv")

"""create Kernel and Kernel matrices"""
for state in np.arange(y.shape[1]):
    currentKernel = Matern52(theta=np.abs(np.random.randn(2)),
                             sigma=np.abs(np.random.randn(1)))
    currentKernel.learnHyperparams(
        currentKernel.theta, #-1000*np.ones(currentKernel.theta.size),
        currentKernel.sigma, # 1,
        y[:, state],
        time,
        normalize=True,
        standardize=standardize,
        newNugget=1e-4,
        anneal=False,
        basinIter=kernelIter
        )
    if not os.path.exists("./hyperparams"):
        os.makedirs("./hyperparams")
    np.savetxt("hyperparams/theta{}.csv".format(state),
               currentKernel.theta)
    np.savetxt("hyperparams/sigma{}.csv".format(state),
               np.asarray(currentKernel.sigma).reshape([1, 1]))
