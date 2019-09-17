# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Class implementing the MCMC sampling scheme used for FGPGM.
"""
import numpy as np
from scipy.stats import norm

def MCMCWithBounds(logTarget, xInit, proposalStds, lowerBounds, upperBounds,
                   nSamples, nBurnin):
    """
    Creates MCMC samples respecting boundaries provided by the experiment.
    Uses zero mean Gaussian distributions as proposal.
    
    Parameters
    ----------
    logTarget:      function taking one vector of length nDim as inputs
                    logarithm of the unnormalized density that should be
                    inferred using the MCMC scheme
    xInit:          vector of length nDim
                    initial values of the chain
    proposalStd:    vector of length nDim
                    Standard deviations for the proposal distributions.
    lowerBounds:    vector of length nDim
                    lower bounds for the arguments of targetDensity
    upperBounds:    vector of length nDim
                    upper bounds for the arguments of targetDensity
    nSamples:       integer
                    amount of MCMC steps to be performed in total
    nBurnin:        integer, smaller than nSamples
                    amount of MCMC steps that should be discarded at the
                    beginning for burn-in.

    Returns
    ----------
    samples:    vector of length nSamples-nBurnin
                samples of the chain, describing the probability density
    """
    xInit = np.asarray(xInit)
    allSamples = []
    xNew = xInit.copy()
    logPNew = logTarget(xInit)
    allSamples.append(xNew)
    nAccepted=np.zeros_like(xInit)
    nRejected=np.zeros_like(xInit)
    for i in np.arange(nSamples+nBurnin):
        currSamples = []
        for dim in range(xInit.size):
            # save state and probability of previous iteration
            xOld = xNew.copy()
            logPOld = logPNew
            # create new state
            proposal = getProposal(mean=xOld[dim],
                                   std=proposalStds[dim],
                                   lowerBound=lowerBounds[dim],
                                   upperBound=upperBounds[dim])
            xPotential = xOld.copy()
            xPotential[dim] = proposal
            # get normalization constantes for proposal distribution
            normalizationOldGivenNew = getValidProbability(
                mean=proposal,
                std=proposalStds[dim],
                lowerBound=lowerBounds[dim],
                upperBound=upperBounds[dim]
                )
            normalizationNewGivenOld = getValidProbability(
                mean=xOld[dim],
                std=proposalStds[dim],
                lowerBound=lowerBounds[dim],
                upperBound=upperBounds[dim]
                )
            # get new function value
            logPPotential = logTarget(xPotential)
            # accept if acceptable
            if logPPotential-logPOld > np.log(normalizationOldGivenNew /
                                              normalizationNewGivenOld):
                # catch potential overflows
                acceptanceProb = 1
            else:
                acceptanceProb = np.exp(logPPotential-logPOld)
                acceptanceProb *= normalizationNewGivenOld \
                                / normalizationOldGivenNew
            if np.random.rand(1) <= acceptanceProb:
                xNew = xPotential.copy()
                if proposal > upperBounds[dim] or proposal < lowerBounds[dim]:
                    pass                    
                    print("Bounds violating value accepted (!)")
                logPNew = logPPotential
                nAccepted[dim] += 1
            else:
                logPNew = logPOld
                nRejected[dim] += 1
                pass
            currSamples.append(xNew)
        if i%1000==0:
            print("iteration {}: {}".format(i, xNew))
        allSamples.append(np.mean(currSamples, axis=0))
    
    return np.asarray(allSamples[nBurnin:]), nAccepted, nRejected

def getProposal(mean, std, lowerBound, upperBound, maxTry=100):
    """
    Samples the Gaussian proposal distribution with given mean and std until
    a sample within the lowerBound and upperBound is found.
    Will throw an error if this takes longer than maxTry tries.
    """
    for i in range(maxTry):
        newX = mean + std*np.random.randn(1)
        if newX < upperBound and newX > lowerBound:
            return np.squeeze(newX)
    raise ValueError("Cannot find proper proposal. Maybe bounds are " + 
        "improperly chosen?")

def getValidProbability(mean, std, lowerBound, upperBound):
    """
    Calculates the probability mass of the Gaussian proposal distribution
    with given mean and std that is within the lower and upper bound.
    """
    return norm.cdf(upperBound, loc=mean, scale=std) \
         - norm.cdf(lowerBound, loc=mean, scale=std)
