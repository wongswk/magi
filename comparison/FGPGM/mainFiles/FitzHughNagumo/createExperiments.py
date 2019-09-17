# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the files used later by FGPGM.
Must be run from the testRun{}.format{i} directories.
"""

import numpy as np
from FGPGM.Experiments.FitzHughNagumo import FHN

from matplotlib import pyplot as plt
plt.switch_backend('agg')

"""seed random number generator and save the seed for reproducability"""
seed = np.random.randint(0, 2**32-1)
np.random.seed(seed)
np.savetxt("seed.csv", np.asarray(seed).reshape([1, 1]))

"""experiment and observation creation"""
# parameters copied from MacDonald overview paper
XInit = [-1, 1]
XInit = np.asarray(XInit)
np.savetxt("XInit.csv", XInit)

theta = [0.2, 0.2, 3]
theta = np.asarray(theta)
np.savetxt("trueODEParams.csv", theta)

print("normal time")
tEnd = 10 # 5
#dt = 0.1 # 0.125
dt = np.loadtxt("dt.csv")
time = np.arange(0, tEnd+0.5*dt, dt)  # obs time
time = np.asarray(time)
np.savetxt("time.csv", time)

SNR = np.loadtxt("SNRInput.csv") # same as VGM and AGM
np.savetxt("usedSNR.csv", np.asarray(SNR).reshape([1, 1]))

standardize=1
np.savetxt("standardize.csv", np.asarray(standardize).reshape([1, 1]))

experiment = FHN()
# sample noisy and true states at observation times
x, y = experiment.sampleTrajectoryNonUniform(XInit, theta, time, obsNoiseStd=None,
                                             SNR=SNR, plotting="./experimentPlots")

np.savetxt("trueStates.csv", x)
np.savetxt("observations.csv", y)
