import numpy as np
from arma import ode_system, solve_gpds


def fOde(theta, x):
    V = x[:, 0]
    R = x[:, 1]
    Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
    Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
    result = np.stack([Vdt, Rdt], axis=1)
    return result

def fOdeDx(theta, x):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDx[:,0,0] = theta[2] * (1 - np.square(V))
    resultDx[:,1,0] = theta[2]
    resultDx[:,0,1] = -1.0 / theta[2]
    resultDx[:,1,1] = -1.0*theta[1]/theta[2]
    return resultDx

def fOdeDtheta(theta, x):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
    resultDtheta[:,0,1] = 1.0 / theta[2]
    resultDtheta[:,1,1] = -R / theta[2]
    resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)
    return resultDtheta

fn_system = ode_system("FN-python", fOde, fOdeDx, fOdeDtheta,
                       thetaLowerBound=np.array([0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf]))

ydataTruth = [[-1, -0.20925154659207, 1.83568745612659,
               2.0135309104164, 1.90873019346171, 1.79486115820876, 1.67211213239152,
               1.53662974499158, 1.38120073986218, 1.19013170106295, 0.919484435213506,
               0.383547960005714, -1.29912297158911, -1.94703008843896, -1.82621515345628,
               -1.67605296462337, -1.50179855978052, -1.28453002771101, -0.968149651951442,
               -0.289583787997096, 1.69708215851227, 2.00615161512914, 1.90282235165016,
               1.78857976062948, 1.6652708465063, 1.52895166419629, 1.37214013119417,
               1.17838132678723, 0.900810540853709, 0.335888855600565, -1.40628542040035,
               -1.94307261810271, -1.81869966439223, -1.66749000485501, -1.49157844936251,
               -1.27103410403028, -0.945810196662687, -0.225460500531252, 1.77161890315082,
               2.00150473412179, 1.89693950309839],
              [1, 1.10971573208142,
               0.973971505460551, 0.645805849456405, 0.335798626053281, 0.0539770521943058,
               -0.199245124193061, -0.423064123771236, -0.615830984225796, -0.774184003197497,
               -0.890478676998256, -0.941793805684236, -0.824808213741384, -0.468260373942925,
               -0.110108123580818, 0.213525581109257, 0.50010099996594, 0.745568908212889,
               0.940528043027411, 1.0548657257884, 0.949537453201733, 0.628529925666569,
               0.320093849068738, 0.0397833888090213, -0.211900330950493, -0.43412034922384,
               -0.625162798467007, -0.78153058734626, -0.895191437797031, -0.941542148713408,
               -0.809587454312223, -0.448430208702261, -0.0919986870657319,
               0.229726996465537, 0.514243895516946, 0.757337477682966, 0.949110759701012,
               1.05705510574894, 0.934260823997242, 0.611444769682256, 0.304475451953246
               ]]
ydataTruth = np.array(ydataTruth).transpose()

ydataV = [-0.86, -0.26, 2.14, 1.94, 1.63, 1.75,
          1.92, 1.39, 1.29, 1.59, 0.63, 0.78, -1.59, -1.92, -1.56, -1.58,
          -1.26, -1.34, -0.62, -0.39, 1.58, 2.29, 1.69, 1.61, 1.88, 1.57,
          1.28, 1.09, 1.21, 0.1, -1.66, -2.05, -1.55, -1.81, -1.72, -0.98,
          -0.77, -0.09, 1.87, 2.18, 1.67]
ydataR = [0.94, 0.87, 0.62, 0.44,
          0.07, 0.02, -0.55, -0.09, -0.66, -0.73, -0.73, -0.63, -0.85,
          -0.55, 0.01, 0.43, 0.4, 0.57, 0.64, 1.26, 1.09, 0.46, 0.13, 0.14,
          -0.3, -0.53, -0.5, -0.35, -1.03, -1.02, -0.6, -0.61, -0.05, 0.31,
          0.82, 0.85, 0.64, 1.31, 0.78, 0.47, 0.35]
ydata = np.stack([np.array(ydataV), np.array(ydataR)], axis=1)
tvecObs = np.linspace(0, 20, num=41)
tvecFull = np.linspace(0, 20, num=161)
yFull = np.ndarray([161, 2])
yFull.fill(np.nan)
yFull[np.linspace(0, 160, num=41).astype(int), :] = ydata

xInitExogenous = np.zeros_like(yFull)
for j in range(2):
    xInitExogenous[:, j] = np.interp(tvecFull, tvecObs, ydata[:, j])

result = solve_gpds(
    yFull,
    fn_system,
    tvecFull,
    sigmaExogenous = np.array([]),
    phiExogenous = np.array([[]]),
    xInitExogenous = xInitExogenous,
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.shape[0]/ydata.shape[0],
    priorTemperatureDeriv = yFull.shape[0]/ydata.shape[0],
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 100,
    burninRatioHmc = 0.5,
    niterHmc = 20001,
    stepSizeFactorHmc = 0.06,
    nEpoch = 1,
    bandSize = 20,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = False,
    verbose = True)

print(result['phiUsed'])
print(result['samplesCpp'])

# verify trajectory RMSE
samplesCpp = result['samplesCpp']
llikId = 0
xId = range(np.max(llikId)+1, np.max(llikId)+yFull.size+1)
thetaId = range(np.max(xId)+1, np.max(xId)+3+1)
sigmaId = range(np.max(thetaId)+1, np.max(thetaId)+yFull.shape[1]+1)

burnin = int(20001*0.5)
xsampled = samplesCpp[xId, (burnin+1):]
xsampled = xsampled.reshape([yFull.shape[1], yFull.shape[0], -1])

from matplotlib import pyplot as plt
for j in range(yFull.shape[1]):
    plt.plot(tvecFull, xsampled[j, :, -1])  # one single sample plot

inferred_trajectory = np.mean(xsampled, axis=-1)
for j in range(yFull.shape[1]):
    plt.plot(tvecFull, inferred_trajectory[j, :])  # inferred trajectory plot


inferred_trajectory_orig_time = np.zeros_like(ydataTruth)
for j in range(yFull.shape[1]):
    inferred_trajectory_orig_time[:, j] = np.interp(tvecObs, tvecFull, inferred_trajectory[j, :])
trajectory_rmse = np.sqrt(np.mean((inferred_trajectory_orig_time - ydataTruth)**2, axis=0))
np.savetxt("trajectory_rmse.csv", trajectory_rmse)
