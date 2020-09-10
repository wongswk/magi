import numpy as np
from pymagi import ArmaVector, ArmaMatrix, ArmaCube, OdeSystem, solveMagiPy


def vector(arma_vector):
    return np.ndarray(buffer=arma_vector, shape=(arma_vector.size(),))


def matrix(arma_matrix):
    return np.ndarray(buffer=arma_matrix, shape=(arma_matrix.n_rows, arma_matrix.n_cols), order="F")


def cube(arma_cube):
    return np.ndarray(buffer=arma_cube, shape=(arma_cube.n_rows, arma_cube.n_cols, arma_cube.n_slices), order="F")


def integer_vector(arma_integer_vector):
    return np.ndarray(buffer=arma_integer_vector, shape=(arma_integer_vector.size(),), dtype=int)


def integer_matrix(arma_integer_matrix):
    return np.ndarray(
        buffer=arma_integer_matrix, shape=(arma_integer_matrix.n_rows, arma_integer_matrix.n_cols), order="F", dtype=int
    )


def ode_system(name, fOde, fOdeDx, fOdeDtheta, thetaLowerBound, thetaUpperBound):
    system = OdeSystem()
    def fOdeArma(theta, x):
        theta = vector(theta)
        x = matrix(x)
        result = fOde(theta, x)
        return ArmaMatrix(result.T.copy())

    def fOdeDxArma(theta, x):
        theta = vector(theta)
        x = matrix(x)
        resultDx = fOdeDx(theta, x)
        return ArmaCube(resultDx.T.copy())

    def fOdeDthetaArma(theta, x):
        theta = vector(theta)
        x = matrix(x)
        resultDtheta = fOdeDtheta(theta, x)
        return ArmaCube(resultDtheta.T.copy())

    system.fOde = fOdeArma
    system.fOdeDx = fOdeDxArma
    system.fOdeDtheta = fOdeDthetaArma
    system.thetaLowerBound = ArmaVector(thetaLowerBound)
    system.thetaUpperBound = ArmaVector(thetaUpperBound)
    system.name = name
    system.thetaSize = thetaLowerBound.size
    return system


def solve_magi(
        yFull,
        odeModel,
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
        verbose = True):

    sigmaExogenous = ArmaVector(np.ndarray(0)) if sigmaExogenous.size == 0 else ArmaVector(sigmaExogenous)
    phiExogenous = ArmaMatrix(np.ndarray([0, 0])) if phiExogenous.size == 0 else ArmaMatrix(phiExogenous).t()
    xInitExogenous = ArmaMatrix(np.ndarray([0, 0])) if xInitExogenous.size == 0 else ArmaMatrix(xInitExogenous).t()
    thetaInitExogenous = ArmaVector(np.ndarray(0)) if thetaInitExogenous.size == 0 else ArmaVector(thetaInitExogenous)
    muExogenous = ArmaMatrix(np.ndarray([0, 0])) if muExogenous.size == 0 else ArmaMatrix(muExogenous).t()
    dotmuExogenous=ArmaMatrix(np.ndarray([0, 0])) if dotmuExogenous.size == 0 else ArmaMatrix(dotmuExogenous).t()

    result_solved = solveMagiPy(
        yFull=ArmaMatrix(yFull).t(),
        odeModel=odeModel,
        tvecFull=ArmaVector(tvecFull),
        sigmaExogenous=sigmaExogenous,
        phiExogenous=phiExogenous,
        xInitExogenous=xInitExogenous,
        thetaInitExogenous=thetaInitExogenous,
        muExogenous=muExogenous,
        dotmuExogenous=dotmuExogenous,
        priorTemperatureLevel=priorTemperatureLevel,
        priorTemperatureDeriv=priorTemperatureDeriv,
        priorTemperatureObs=priorTemperatureObs,
        kernel=kernel,
        nstepsHmc=nstepsHmc,
        burninRatioHmc=burninRatioHmc,
        niterHmc=niterHmc,
        stepSizeFactorHmc=stepSizeFactorHmc,
        nEpoch=nEpoch,
        bandSize=bandSize,
        useFrequencyBasedPrior=useFrequencyBasedPrior,
        useBand=useBand,
        useMean=useMean,
        useScalerSigma=useScalerSigma,
        useFixedSigma=useFixedSigma,
        verbose=verbose)

    phiUsed = matrix(result_solved.phiAllDimensions)
    phiUsed = np.copy(phiUsed.reshape([-1])).reshape([2, -1])
    samplesCpp = matrix(result_solved.llikxthetasigmaSamples.slice(0))
    return dict(phiUsed=phiUsed, samplesCpp=samplesCpp)
