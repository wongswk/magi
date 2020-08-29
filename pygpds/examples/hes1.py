import numpy as np
from arma import ode_system, solve_gpds


def fOde(theta, x):
    P = np.exp(x[:, 0])
    M = np.exp(x[:, 1])
    H = np.exp(x[:, 2])
    PMHdt = np.zeros(shape=np.shape(x))
    PMHdt[:,0] = -theta[0]*H + theta[1]*M/P - theta[2]
    PMHdt[:,1] = -theta[3] + theta[4]/(1+np.square(P))/M
    PMHdt[:,2] = -theta[0]*P + theta[5]/(1+np.square(P))/H - theta[6]
    return PMHdt


def fOdeDx(theta, x):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    P = (x[:, 0])
    M = (x[:, 1])
    H = (x[:, 2])

    expMminusP = np.exp(M-P)
    dP = -np.power(1+np.exp(2*P), -2) * np.exp(2*P)*2

    resultDx[:,0,0] = -theta[1]*expMminusP
    resultDx[:,1,0] = theta[1]*expMminusP
    resultDx[:,2,0] = -theta[0]*np.exp(H)

    resultDx[:,0,1] = theta[4]*np.multiply(np.exp(-M),dP)
    resultDx[:,1,1] = -theta[4]*np.exp(-M)/(1+np.exp(2*P))

    resultDx[:,0,2] = -theta[0]*np.exp(P) + theta[5]*np.multiply(np.exp(-H),dP)
    resultDx[:,2,2] = -theta[5]*np.exp(-H)/(1+np.exp(2*P))

    return resultDx


def fOdeDtheta(theta, x):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    P = (x[:, 0])
    M = (x[:, 1])
    H = (x[:, 2])

    resultDtheta[:,0,0] = -np.exp(H)
    resultDtheta[:,1,0] = np.exp(M-P)
    resultDtheta[:,2,0].fill(-1)

    resultDtheta[:,3,1].fill(-1)
    resultDtheta[:,4,1] = np.exp(-M)/(1+np.exp(2*P))

    resultDtheta[:,0,2] = -np.exp(P)
    resultDtheta[:,5,2] = np.exp(-H)/(1+np.exp(2*P))
    resultDtheta[:,6,2].fill(-1)

    return resultDtheta


hes1_system = ode_system("Hes1-log-python", fOde, fOdeDx, fOdeDtheta,
                         thetaLowerBound=np.ones(7) * 0, thetaUpperBound=np.ones(7) * np.inf)


ydataP = [0.74, 0.78, 1.86, 1.86, 2.2, 1.93, 1.47, 1.03, 0.36,
          0.88, 1.68, 1.97, 2.15, 1.85, 1.8, 1.47, 0.71]
ydataM = [0.91, 0.82, 0.71, -0.11, 0.08, -0.45, -0.05, 0.2,
          0.88, 1.09, 0.3, 0.35, 0.25, -0.23, -0.51, -0.09]
yFull = np.ndarray([33, 3])
yFull.fill(np.nan)
yFull[np.linspace(0, 32, num=17).astype(int), 0] = ydataP
yFull[np.linspace(1, 31, num=16).astype(int), 1] = ydataM

tvecFull = np.linspace(0, 240, num=33)


result = solve_gpds(
    yFull,
    hes1_system,
    tvecFull,
    sigmaExogenous = np.array([0.15, 0.15, 0.15]),
    phiExogenous = np.array([[]]),
    xInitExogenous = np.array([[]]),
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureDeriv = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 500,
    burninRatioHmc = 0.5,
    niterHmc = 20001,
    stepSizeFactorHmc = 0.01,
    nEpoch = 1,
    bandSize = 20,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = True,
    verbose = True)

print(result['phiUsed'])
print(result['samplesCpp'])
