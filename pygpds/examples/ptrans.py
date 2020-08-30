import numpy as np
from arma import ode_system, solve_gpds


def fOde(theta, x):
    resultdt = np.zeros(shape=np.shape(x))

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultdt[:, 0] = -theta[0]*S - theta[1] * S @ R + theta[2] * RS
    resultdt[:, 1] = theta[0]*S
    resultdt[:, 2] = -theta[1]*S@R + theta[2]*RS + theta[4] * RPP / (theta[5]+RPP)
    resultdt[:, 3] = theta[1]*S@R - theta[2]* RS - theta[3]*RS
    resultdt[:, 4] = theta[3]*RS - theta[4] * RPP / (theta[5]+RPP)

    return resultdt


def fOdeDx(theta, x):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultDx[:, 0, 0] = -theta[0] - theta[1] * R
    resultDx[:, 2, 0] = -theta[1] * S
    resultDx[:, 3, 0].fill(theta[2])

    resultDx[:, 0, 1].fill(theta[0])

    resultDx[:, 0, 2] = -theta[1]*R
    resultDx[:, 2, 2] = -theta[1]*S
    resultDx[:, 3, 2].fill(theta[2])
    resultDx[:, 4, 2] =  theta[4] * theta[5] /  np.square(theta[5] + RPP)

    resultDx[:, 0, 3] = theta[1]*R
    resultDx[:, 2, 3] = theta[1]*S
    resultDx[:, 3, 3].fill(-theta[2] - theta[3])

    resultDx[:, 3, 4].fill(theta[3])
    resultDx[:, 4, 4] = -theta[4] * theta[5] /  np.square(theta[5] + RPP)

    return resultDx


def fOdeDtheta(theta, x):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultDtheta[:, 0, 0] = -S
    resultDtheta[:, 1, 0] = -S@R
    resultDtheta[:, 2, 0] = RS

    resultDtheta[:, 0, 1] = S

    resultDtheta[:, 1, 2] = -S@R
    resultDtheta[:, 2, 2] = RS
    resultDtheta[:, 4, 2] = RPP / (theta[5]+RPP)
    resultDtheta[:, 5, 2] = -theta[4] * RPP / np.square(theta[5]+RPP)

    resultDtheta[:, 1, 3] = S@R
    resultDtheta[:, 2, 3] = -RS
    resultDtheta[:, 3, 3] = -RS

    resultDtheta[:, 3, 4] = RS
    resultDtheta[:, 4, 4] = - RPP / (theta[5]+RPP)
    resultDtheta[:, 5, 4] = theta[4] * RPP / np.square(theta[5]+RPP)

    return resultDtheta


ptrans_system = ode_system("PTrans-python", fOde, fOdeDx, fOdeDtheta,
                           thetaLowerBound=np.ones(6) * 0, thetaUpperBound=np.ones(6) * 4)


ydata = [[1.00003, 0.5887, 0.40509, 0.23316, 0.18567,
          0.12149, 0.06848, 0.024, 0.00674, 0.00109, 0.00044, -0.00099,
          -0.00039, 0.00012, 0.00061],
         [0.00046, 0.05223, 0.08721,
          0.1315, 0.14517, 0.16616, 0.18425, 0.19958, 0.20584, 0.20627,
          0.2076, 0.20519, 0.20674, 0.20635, 0.20808],
         [0.99956,
          0.64249, 0.49785, 0.38488, 0.36059, 0.3376, 0.33559, 0.36219,
          0.40655, 0.51279, 0.61291, 0.70233, 0.77953, 0.89698, 0.95844
          ],
         [0.00244, 0.3011, 0.35024, 0.2845, 0.23995, 0.16463,
          0.08887, 0.03159, 0.01217, -0.00032, -0.00066, -0.00087, 0.00052,
          0.00104, 0.00054],
         [-0.00172, 0.05523, 0.15238, 0.33098,
          0.40037, 0.49987, 0.5774, 0.60517, 0.57964, 0.48722, 0.38806,
          0.29898, 0.21757, 0.1039, 0.04104]]
ydata = np.array(ydata).transpose()
tvecObs = [0, 1, 2, 4, 5, 7, 10, 15, 20, 30, 40, 50, 60, 80, 100]
yFillI0 = np.ndarray([101, 5])
yFillI0.fill(np.nan)
yFillI0[tvecObs, :] = ydata
tvecI0 = np.linspace(0, 100, num=101)
for j in range(5):
    yFillI0[:, j] = np.interp(tvecI0, tvecObs, ydata[:, j])


yFull = np.ndarray([201, 5])
yFull.fill(np.nan)
yFull[tvecObs, :] = ydata
tvecFull = np.linspace(0, 200, num=201)
xInitExogenous = np.zeros_like(yFull)
for j in range(5):
    xInitExogenous[:, j] = np.interp(tvecFull, tvecObs, ydata[:, j])


hyperInit = solve_gpds(
    yFillI0,
    ptrans_system,
    tvecI0,
    sigmaExogenous = np.array([]),
    phiExogenous = np.array([[]]),
    xInitExogenous = np.array([[]]),
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureDeriv = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 100,
    burninRatioHmc = 0.5,
    niterHmc = 2,
    stepSizeFactorHmc = 0.01,
    nEpoch = 1,
    bandSize = 40,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = False,
    verbose = True)

print(hyperInit['phiUsed'])
print(hyperInit['samplesCpp'])
sigmaUsed = hyperInit['samplesCpp'][-5:, 0]

# NEW (Aug 11) ----- plug in sigma estimate and re-estimate phi
hyperInit = solve_gpds(
    yFillI0,
    ptrans_system,
    tvecI0,
    sigmaExogenous = sigmaUsed,
    phiExogenous = np.array([[]]),
    xInitExogenous = np.array([[]]),
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureDeriv = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 100,
    burninRatioHmc = 0.5,
    niterHmc = 2,
    stepSizeFactorHmc = 0.01,
    nEpoch = 1,
    bandSize = 40,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = False,
    verbose = True)
phiUsed = hyperInit['phiUsed']

# sampling
result = solve_gpds(
    yFull,
    ptrans_system,
    tvecFull,
    sigmaExogenous = sigmaUsed,
    phiExogenous = phiUsed,
    xInitExogenous = xInitExogenous,
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureDeriv = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 100,
    burninRatioHmc = 0.5,
    niterHmc = 20001,
    stepSizeFactorHmc = 0.01,
    nEpoch = 1,
    bandSize = 40,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = False,
    verbose = True)

print(result['phiUsed'])
print(result['samplesCpp'])
