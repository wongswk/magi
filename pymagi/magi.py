import numpy as np
from arma import ode_system, solve_magi

def MagiSolver(y, odeModel, tvec=None, control=dict()):
    if tvec is None:
        tvec = y[:, 0]
        y = y[:, 1:]
        print('tvec is not specified, first column in y will be used as time')

    if 'sigma' in control.keys():
        sigmaExogenous = control['sigma']
        sigmaExogenous[~np.isfinite(sigmaExogenous)] = 1
    else:
        sigmaExogenous = np.array([])

    if 'phi' in control.keys():
        phiExogenous = control['phi']
    else:
        phiExogenous = np.array([])

    if 'xInit' in control.keys():
        xInitExogenous = control['xInit']
    else:
        xInitExogenous = np.array([])

    if 'thetaInit' in control.keys():
        thetaInitExogenous = control['thetaInit']
    else:
        thetaInitExogenous = np.array([])

    if 'mu' in control.keys():
        muExogenous = control['mu']
    else:
        muExogenous = np.array([])

    if 'dotmu' in control.keys():
        dotmuExogenous = control['dotmu']
    else:
        dotmuExogenous = np.array([])

    if 'priorTemperature' in control.keys():
        priorTemperatureLevel = control['priorTemperature']
        priorTemperatureDeriv = control['priorTemperature']
    else:
        priorTemperatureLevel = 1/np.mean(np.isfinite(y))
        priorTemperatureDeriv = 1/np.mean(np.isfinite(y))

    if 'dotmu' in control.keys():
        dotmuExogenous = control['dotmu']
    else:
        dotmuExogenous = np.array([])

    if 'niterHmc' in control.keys():
        niterHmc = control['niterHmc']
    else:
        niterHmc = 20000

    if 'burninRatio' in control.keys():
        burninRatio = control['burninRatio']
    else:
        burninRatio = 0.5

    if 'nstepsHmc' in control.keys():
        nstepsHmc = control['nstepsHmc']
    else:
        nstepsHmc = 200

    if 'stepSizeFactor' in control.keys():
        stepSizeFactor = control['stepSizeFactor']
    else:
        stepSizeFactor = 0.01

    if 'bandSize' in control.keys():
        bandSize = control['bandSize']
    else:
        bandSize = 20

    if 'useFixedSigma' in control.keys():
        useFixedSigma = control['useFixedSigma']
    else:
        useFixedSigma = False


    result = solve_magi(
        y,
        odeModel,
        tvec,
        sigmaExogenous = sigmaExogenous,
        phiExogenous = phiExogenous,
        xInitExogenous = xInitExogenous,
        thetaInitExogenous = thetaInitExogenous,
        muExogenous = muExogenous,
        dotmuExogenous = dotmuExogenous,
        priorTemperatureLevel = priorTemperatureLevel,
        priorTemperatureDeriv = priorTemperatureDeriv,
        priorTemperatureObs = 1.0,
        kernel = "generalMatern",
        nstepsHmc = nstepsHmc,
        burninRatioHmc = burninRatio,
        niterHmc = niterHmc,
        stepSizeFactorHmc = stepSizeFactor,
        nEpoch = 1,
        bandSize = bandSize,
        useFrequencyBasedPrior = True,
        useBand = True,
        useMean = True,
        useScalerSigma = False,
        useFixedSigma = useFixedSigma,
        verbose = True)

    phiUsed = result['phiUsed']
    samplesCpp = result['samplesCpp']

    llikId = 0
    xId = range(np.max(llikId)+1, np.max(llikId)+y.size+1)
    thetaId = range(np.max(xId)+1, np.max(xId)+3+1)
    sigmaId = range(np.max(thetaId)+1, np.max(thetaId)+yFull.shape[1]+1)

    burnin = int(niterHmc*burninRatio)
    xsampled = samplesCpp[xId, (burnin+1):]
    xsampled = xsampled.reshape([y.shape[1], y.shape[0], -1])

    thetaSampled = samplesCpp[thetaId, (burnin+1):]

    return dict(
        theta=thetaSampled,
        xsampled=xsampled,
        lp=samplesCpp[llikId, (burnin+1):],
        sigma=samplesCpp[sigmaId, (burnin+1):],
        phi=phiUsed,
    )
