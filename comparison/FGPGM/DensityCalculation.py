# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Helper functions to calculate matrices and densities as described in the
FGPGM paper.
"""
import numpy as np

def getAs(CDashs, DashCs, CPhis, CDoubleDashs):
    """
    Parameters
    ----------
    CDashs:         list of matrices of shape nTime x nTime
                    each entry represents a state. matrices as returned by
                    kernel.getC_PhiDash
    DashCs:         list of matrices of shape nTime x nTime
    CInvs:          list of matrices of shape nTime x nTime
    CDoubleDashs:   list of matrices of shape nTime x nTime

    Returns
    ----------
    A:  list of matrices of shape nTime x nTime
        each entry represents one state
    """
    A = []
    for i in np.arange(len(CDashs)):
        A.append(
            CDoubleDashs[i] - np.dot(
                DashCs[i],
                np.linalg.solve(CPhis[i], CDashs[i])))
    return A

def getDs(DashCs, CPhis):
    """
    each entry represents a state
    
    Parameters
    ----------
    DashCs:         list of matrices of shape nTime x nTime
    CInvs:          list of matrices of shape nTime x nTime

    Returns
    ----------
    D:  list of matrices of shape nTime x nTime
        each entry represents one state
    """
    D = []
    for i in np.arange(len(DashCs)):
        def getProdWithD(x, index=i):
            """
            defines a function to get a product of a vector with the matrix D
            """
            return np.dot(DashCs[index],
                          np.linalg.solve(CPhis[index], x)
                          )
        D.append(getProdWithD)
    return D

def getLambdaStars(gamma, nTime):
    """
    each entry represents the LambdaStar of one state, which is the noise
    covariance between the ODEs and the GP estimates of the derivatives
    
    Parameters
    ----------
    gamma:  vector of length nStates
            noise estimate on the derivatives
    nTime:  scalar
            amount of time steps in this experiment
    Returns
    ----------
    Lambdas:    list of nStates matrices of shape nTime x nTime    
    """
    gamma = np.asarray(gamma)
    Lambdas = []
    for i in np.arange(gamma.shape[0]):
        Lambdas.append(np.eye(nTime)*gamma[i])
    return Lambdas

def calculateGPPostLog(yNormal, xPret, CPhis, sigmas, mean=None, std=None,
                       includeDet=False):
    """
    calculates the logarithm of the GP posterior of the states
    
    y: unfolded observations
    x: unfolded states, pretreated according to normalization/standardization
    mean: unfolded means
    std:  unfolded stds
    CPhis: list of covariance matrices. Stacked results of kernel.getCPhi(time)
    sigmas: vector of length nStates with the (estimated) stds of the
            observation noise
    """
    yNormal = np.asarray(yNormal)
    xPret = np.asarray(xPret)
    sigmas = np.asarray(sigmas)
    assert xPret.shape[0] == yNormal.shape[0]
    assert len(CPhis) == sigmas.shape[0]
    nTime = CPhis[0].shape[0]
    contributions = np.zeros(len(CPhis))

    yPret = (yNormal - mean) / std

    # iterate through states and calculate respective GP posterior equivalents        
    for i in np.arange(len(CPhis)):
        # iterate through all states
        startIndex = nTime*i
        endIndex = nTime*(i+1)
        currentXPret = xPret[startIndex:endIndex]
        currentYPret = yPret[startIndex:endIndex]

        priorContrib = -1./2*np.dot(
            currentXPret,
            np.linalg.solve(CPhis[i], currentXPret)
            )
        if includeDet:
            detArgumentPrior = -1./2*np.prod(np.linalg.slogdet(CPhis[i]))
            priorContrib += detArgumentPrior

        difference = currentXPret - currentYPret

        obsContrib = -sigmas[i]**(-2)/2*np.dot(difference, difference)
        if includeDet:
            detArgumentObs = -1./2*np.prod(np.linalg.slogdet(
                sigmas[i]**2*np.eye(CPhis[i].shape[0])))
            obsContrib += detArgumentObs
        contributions[i] = priorContrib + obsContrib
    return np.sum(contributions)

def calculateFPostLog(f, xPret, theta, As, Ds, Lambdas, mean=None, std=None,
                      includeDet=False):
    """
    calculates the second component of the density function, which corresponds
    to the logarithm of the posterior on F.

    Parameters
    ----------
    f:  function handle
        function representing the ODEs by mapping x[t] and theta[t] to x_dot[t]
        takes x as first and theta as second argument
    xPret:  vector of shape nTime*nStates
        unfolded states. [x1[t0], x1[t1], ..., x1[tEnd], x2[t0]...]
        pretreated as specified by normalize/standardize
    mean:   Vector of same shape as x
            if None, it is assumed that no GP normalization has been done
            if not None, it is expected to contain the mean of each batch of
            observation as entries.
    std:    Vector of same shape as x
            if None, it is assumed that no standardization has been done
            if not None, it is assumed to contain the std of each batch of
            observations as entries
    """
    xPret = np.asarray(xPret)
    assert len(As) == len(Ds)
    assert len(Ds) == len(Lambdas)
    nTime = As[0].shape[0]
    nStates = len(As)
    assert nStates*nTime == xPret.shape[0]

    # calculate true states and corrected unfolded one
    mean = np.asarray(mean)
    std = np.asarray(std)
    assert mean.size == xPret.size
    assert std.size == xPret.size
    xNormal = std*xPret + mean

    # refold the states to get a matrix of shape nStates x nTime
    xMatrixPret = np.zeros([nStates, nTime])
    xMatrixNormal = np.zeros([nStates, nTime])
    stdMatrix = np.zeros([nStates, nTime])
    for i in np.arange(nStates):
        xMatrixPret[i, :] = xPret[i*nTime:(i+1)*nTime]
        xMatrixNormal[i, :] = xNormal[i*nTime:(i+1)*nTime]
        stdMatrix[i, :] = std[i*nTime:(i+1)*nTime]

    # calculate derivatives
    fMatrix = np.zeros_like(xMatrixPret)
    probabilities = np.zeros(nStates)

    for t in np.arange(nTime):
        fMatrix[:, t] = f(xMatrixNormal[:, t], theta)

    # normalize derivatives
    fMatrix = fMatrix / stdMatrix

    for state in np.arange(nStates):
        currentMean = Ds[state](xMatrixPret[state, :])
        currentDiff = fMatrix[state, :] - currentMean
        currentSigma = As[state] + Lambdas[state]

        probabilities[state] = -1./2*np.dot(
            currentDiff,
            np.linalg.solve(currentSigma, currentDiff)
            )
        if includeDet:
            detArgumentProb = -1./2*np.prod(np.linalg.slogdet(currentSigma))
            probabilities[i] += detArgumentProb

    return np.sum(probabilities)

def calculateLogDensity(y, x, CPhis, sigmas, f, theta, As, Ds, Lambdas,
                        mean=None, std=None, includeDet=False):
    """
    calculates a function proportional to the log density for given inputs
    
    normalize:  array of shape like y or None
                means of the observations used in calculating the hyperparams.
                if None, it is assumed that no normalization has been done
    standardize: array of shape like y or None
                stds of the observations used in calculating the hyperparams.
                if None, it is assumed that no normalization has been done

    y and x are the observations and the states unfolded like
    [x1[t0], x1[t1], ..., x1[tEnd], x2[t0]...]
    if either normalize or standardize, x is assumed to ALREADY BE TRANSFORMED
    in that way
    
    f is the function representing the odes, which takes x as first and
    theta as second input.
    
    theta are the parameters for the ODE
    
    sigmas is a vector of length n containing the observation noise per state

    CPhis is a list of length n containing the covariance matrices of the
    kernels
    """
    scale = 2*float(x.size)
    y = np.asarray(y)
    assert len(As) == len(Ds)
    assert len(Ds) == len(Lambdas)

    if mean is None:
        mean = np.zeros_like(y)
    if std is None:
        std = np.ones_like(y)

    prob = calculateFPostLog(
                             f, x, theta, As, Ds, Lambdas, mean=mean, std=std,
                             includeDet=includeDet) \
         + calculateGPPostLog(y, x, CPhis, sigmas, mean=mean, std=std,
                              includeDet=includeDet)
    return -prob / scale

def _getMultivariateDensityLog(x, mu, Sigma=None, Gamma=None, includeDet=False):
    """
    Calculates the important parts of the logarithm of a multivariate Gaussian
    density.
    
    Parameters
    ----------
    x:          vector of length nStates
                point at which log of multivariate density should be evaluated
    mu:         vector of length nStates
                mean of the multivariate density
    Sigma:      matrix of shape nStates x nStates or None
                covariance matrix of the pdf. Can only be None if Gamma is not
                None. Will be ignored if Gamma is not None.
    Gamma:      matrix of shape nStates x nStates or None
                precision matrix of the pdf. Can only be None if Sigma is not
                None. If Gamma is not None, Sigma will be ignored
    includeDet: boolean
                if true, the determinant of the log of the pdf will be returned
                as well.
                if false, only the exponent of the exponential will be returned.
                The constant term is always ignored.
    """
    assert ((not Gamma is None) or (not Sigma is None))
    x = np.asarray(x)
    difference = x - mu

    if Gamma is not None:
        expArgument = -1./2*np.linalg.multi_dot([difference, Gamma, difference])
    else:
        expArgument = np.linalg.solve(Sigma, difference)
        expArgument = -1./2*np.dot(difference, expArgument)

    if includeDet:
        if Gamma is not None:
            detArgument = 1./2*np.prod(np.linalg.slogdet(Gamma))
        else:
            detArgument = -1./2*np.prod(np.linalg.slogdet(Sigma))
        return expArgument + detArgument
    else:
        return expArgument