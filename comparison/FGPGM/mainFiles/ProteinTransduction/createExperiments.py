# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the files used later.
Must be run from a folder where the output should be saved.
"""
import numpy as np
from FGPGM.Experiments.ProteinTransduction import ProteinTransduction as PT

from matplotlib import pyplot as plt
plt.switch_backend('agg')

"""seed random number generator and save the seed for reproducability"""
seed = np.random.randint(0, 2**32-1)
np.random.seed(seed)
np.savetxt("seed.csv", np.asarray(seed).reshape([1, 1]))

"""experiment and observation creation"""
# initial values
XInit = [1,  # initial values
         0,
         1,
         0,
         0]
XInit = np.asarray(XInit)
np.savetxt("XInit.csv", XInit)
# parameters of the ODEs
theta = [0.07,  # ODE parameters
         0.6,
         0.05,
         0.3,
         0.017,
         .3]
theta = np.asarray(theta)
np.savetxt("trueODEParams.csv", theta)
# observation points of the experiment
print("normal time")
time = [0, 1, 2, 4, 5, 7, 10, 15, 20, 30, 40, 50, 60, 80, 100]  # obs time
time = np.asarray(time)
np.savetxt("time.csv", time)
# standard deviation of the experiment
noiseObsStd=0.01
np.savetxt("obsStdUsed.csv", np.asarray(noiseObsStd).reshape([1, 1]))
# create experiment
experiment = PT()
# sample noisy and true states at observation times
plotPath=None
x, y = experiment.sampleTrajectoryNonUniform(XInit,
                                             theta,
                                             time,
                                             noiseObsStd,
                                             plotting=plotPath
                                             )
np.savetxt("trueStates.csv", x)
np.savetxt("observations.csv", y)