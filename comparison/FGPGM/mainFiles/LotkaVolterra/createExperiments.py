# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the files used later.
Must be run from a folder where the output should be saved.
"""
import numpy as np
from FGPGM.Experiments.LotkaVolterra import LotkaVolterra as LV

from matplotlib import pyplot as plt
plt.switch_backend('agg')

"""seed random number generator and save the seed for reproducability"""
seed = np.random.randint(0, 2**32-1)
np.random.seed(seed)
np.savetxt("seed.csv", np.asarray(seed).reshape([1, 1]))

"""experiment and observation creation - Setup copied from AGM/VGM"""
# initial values
XInit = np.asarray([5, 3])
np.savetxt("XInit.csv", XInit)
# parameters of the ODEs
theta = np.asarray([2, 1, 4, 1])
theta = np.asarray(theta)
np.savetxt("trueODEParams.csv", theta)
# observation points of the experiment
tEnd = 2
dt = 0.1
time = np.arange(0, tEnd+0.5*dt, dt)  # obs time
time = np.asarray(time)
np.savetxt("time.csv", time)
# standard deviation of the experiment
obsStd = 0.1
np.savetxt("obsStdInput.csv", np.asarray(obsStd).reshape([1, 1]))
# create experiment
experiment = LV()
# sample noisy and true states at observation times
plotPath = None
x, y = experiment.sampleTrajectoryNonUniform(XInit,
                                             theta,
                                             time,
                                             obsNoiseStd=obsStd,
                                             SNR=None,
                                             plotting=plotPath)
np.savetxt("trueStates.csv", x)
np.savetxt("observations.csv", y)