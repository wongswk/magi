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
tvecFull = np.linspace(0, 20, num=81)
yFull = np.ndarray([81, 2])
yFull.fill(np.nan)
yFull[np.linspace(0, 80, num=41).astype(int), :] = ydata

result = solve_gpds(
    yFull,
    fn_system,
    tvecFull,
    sigmaExogenous = np.array([]),
    phiExogenous = np.array([[]]),
    xInitExogenous = np.array([[]]),
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = 1.0,
    priorTemperatureDeriv = 1.0,
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 100,
    burninRatioHmc = 0.5,
    niterHmc = 1000,
    stepSizeFactorHmc = 0.1,
    nEpoch = 1,
    bandSize = 20,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = False,
    useScalerSigma = False,
    useFixedSigma = False,
    verbose = True)

print(result['phiUsed'])
print(result['samplesCpp'])
