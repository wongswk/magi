import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pymagi import ArmaVector, ArmaMatrix, ArmaCube, OdeSystem, solveMagiPy, gpsmooth


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
    def fOdeArma(theta, x, tvec):
        theta = vector(theta)
        x = matrix(x)
        tvec = vector(tvec)
        result = fOde(theta, x, tvec)
        return ArmaMatrix(result.T.copy())

    def fOdeDxArma(theta, x, tvec):
        theta = vector(theta)
        x = matrix(x)
        tvec = vector(tvec)
        resultDx = fOdeDx(theta, x, tvec)
        return ArmaCube(resultDx.T.copy())

    def fOdeDthetaArma(theta, x, tvec):
        theta = vector(theta)
        x = matrix(x)
        tvec = vector(tvec)
        resultDtheta = fOdeDtheta(theta, x, tvec)
        return ArmaCube(resultDtheta.T.copy())

    system.fOde = fOdeArma
    system.fOdeDx = fOdeDxArma
    system.fOdeDtheta = fOdeDthetaArma
    system.thetaLowerBound = ArmaVector(thetaLowerBound)
    system.thetaUpperBound = ArmaVector(thetaUpperBound)
    system.name = name
    system.thetaSize = thetaLowerBound.size
    return system

def testDynamicalModel(modelODE, modelDx, modelDtheta, modelName, x, theta, tvec):
    deltaSmall = 1e-06
    tolerance = 1e-04

    f = modelODE(theta, x, tvec)
    fdX = modelDx(theta, x, tvec)

    numericalDx = np.ndarray(shape=[np.shape(tvec)[0], np.shape(x)[1], np.shape(x)[1]]);
    numericalDx.fill(np.nan)
    for j in range(np.shape(x)[1]):
        xnew = x.copy()
        xnew[:,j] = xnew[:,j] + deltaSmall
        numericalDx[:,j,:] = (modelODE(theta, xnew, tvec) - f) / deltaSmall
    
    passedDx = np.amax( np.abs(numericalDx - fdX)) < tolerance
    
    fDtheta = modelDtheta(theta, x, tvec)
    
    numericalDtheta = np.ndarray(shape=[np.shape(tvec)[0], np.shape(theta)[0], np.shape(x)[1]])
    numericalDtheta.fill(np.nan)
    
    for j in range(np.shape(theta)[0]):
        thetanew = theta.copy()
        thetanew[j] = thetanew[j] + deltaSmall
        numericalDtheta[:,j,:] = (modelODE(thetanew, x, tvec) - f) / deltaSmall
    
    passedDtheta = np.amax( np.abs(numericalDtheta - fDtheta))  < tolerance
    
    return(dict(modelName = modelName, Dx = passedDx, Dtheta = passedDtheta))
        
def setDiscretization(dat, level = 0, by = None):
    if by is None:
        if level == 0:
            return dat
        newdat = dat.copy()        
        newdat = newdat[np.argsort(newdat[:,0])]
        dumdat = newdat[1:,:].copy()
        dumdat.fill(np.nan)
        dumdat[:,0] = (newdat[1:,0] + newdat[:-1,0])/2
        newdat = np.append(newdat, dumdat, axis = 0)
        newdat = newdat[np.argsort(newdat[:,0])]
        return setDiscretization(newdat, level = level - 1)
    
    if by > 0:
        tvec = np.arange(np.amin(dat[:,0]), np.amax(dat[:,0]) + by, by)
        tvecAdd = np.setdiff1d(tvec, np.intersect1d(tvec, dat[:,0]))
        newdat = np.ndarray(shape=[np.shape(tvecAdd)[0], np.shape(dat)[1]])
        newdat.fill(np.nan)
        newdat[:,0] = tvecAdd
        newdat = np.append(newdat, dat, axis = 0)
        newdat = newdat[np.argsort(newdat[:,0])]
        return newdat;

def gpsmoothing(yobs, tvec, kerneltype = 'generalMatern', sigma = None):
    distInput = np.abs(tvec[:, None] - tvec)
    yInput = yobs - np.mean(yobs)

    if sigma is None:
        res = gpsmooth(yobsInput = ArmaMatrix(yInput.reshape(1,-1)),
                       distInput = ArmaMatrix(distInput),
                       kernelInput = kerneltype,
                       sigmaExogenScalar = -1.0,
                       useFrequencyBasedPrior = True)
        return dict(sigma = res[2], phi = [res[0], res[1]])

    if sigma > 0:
        res = gpsmooth(yobsInput = ArmaMatrix(yInput.reshape(1,-1)),
                       distInput = ArmaMatrix(distInput),
                       kernelInput = kerneltype,
                       sigmaExogenScalar = sigma,
                       useFrequencyBasedPrior = True)
        return dict(sigma = sigma, phi = [res[0], res[1]])

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

def summaryMagiOutput(x, par_names, sigma = False, lower = 0.025, upper = 0.975):
    
    if sigma:
        allpar = np.vstack((x['theta'], x['sigma']))
    else:
        allpar = x['theta']
    
    theta_est = np.vstack((np.mean(allpar, axis=-1),
                   np.quantile(allpar, q=[lower, upper], axis=-1)))
    return(pd.DataFrame(theta_est, columns=par_names, index=['Mean', f'{100*lower}%', f'{100*upper}%']))

def plotMagiOutput(x, comp_names, obs = True, ci = True, lower = 0.025, upper = 0.975, legend = True):
    xMean = np.mean(x['xsampled'], axis=-1)
    xLB = np.quantile(x['xsampled'], lower, axis=-1)
    xUB = np.quantile(x['xsampled'], upper, axis=-1)
    
    for j in range(xMean.shape[0]):
        plt.plot(x['tvec'], xMean[j,:], label='inferred trajectory')
        
        if ci:
            plt.fill_between(x['tvec'], xLB[j,:], xUB[j,:],
                            alpha=0.2, label=f'{100*(upper-lower)}% credible interval')
        if obs:
            plt.scatter(x['tvec'], x['y'][:,j], label='noisy observations', c='black')
        plt.title(comp_names[j])
        plt.xlabel('Time')
        
        if legend:
            plt.legend()
        plt.show()

