# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Creates and saves all the files used later by FGPGM.
Must be run from the testRun{}.format{i} directories.
"""

import argparse
import numpy as np
from FGPGM.Experiments.FitzHughNagumo import FHN

from matplotlib import pyplot as plt
plt.switch_backend('agg')


"""parse the command line argument"""
parser = argparse.ArgumentParser()
parser.add_argument('--slurm_array_task_id', help='slurm_array_task_id of job array on odyssey', default=-1, type=int)
args = parser.parse_args()


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
dt_candidates = [0.5, 0.4, 0.25, 0.2, 0.1]
if args.slurm_array_task_id < 0:
    dt = 0.25 # 0.125
else:
    dt = dt_candidates[args.slurm_array_task_id % len(dt_candidates)]
    args.slurm_array_task_id = args.slurm_array_task_id // len(dt_candidates)
time = np.arange(0, tEnd+0.5*dt, dt)  # obs time
time = np.asarray(time)
np.savetxt("time.csv", time)

standardize=1
np.savetxt("standardize.csv", np.asarray(standardize).reshape([1, 1]))

# obsStd = 0.2
# np.savetxt("obsStdInput.csv", np.asarray(obsStd).reshape([1, 1]))

# SNR = np.loadtxt("SNRInput.csv")
snr_candidates = [10, 100]
if args.slurm_array_task_id < 0:
    SNR = 10
else:
    SNR = snr_candidates[args.slurm_array_task_id % len(snr_candidates)]
np.savetxt("usedSNR.csv", np.asarray(SNR).reshape([1, 1]))

experiment = FHN()
# sample noisy and true states at observation times
x, y = experiment.sampleTrajectoryNonUniform(XInit, theta, time, obsNoiseStd=None,
                                             SNR=SNR, plotting=None)

np.savetxt("trueStates.csv", x)
np.savetxt("observations.csv", y)