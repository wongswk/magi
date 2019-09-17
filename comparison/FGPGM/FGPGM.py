# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Main class implementing the FGPGM algorithm.
"""
import numpy as np

from .DensityCalculation import getAs, getDs, getLambdaStars, \
     calculateLogDensity

from .MCMCSampler import MCMCWithBounds


class FGPGM(object):
    
    def __init__(self, kernels, time, y, experiment, nODEParams, gamma=1e-3,
                 normalize=False, standardize=False):
        """
        Parameters
        ----------
        kernels:    list of FGPGM.Kernel objects
                    kernels with pretrained hyperparameters
        time:       vector of length nTime
                    time points at which the observations have been made
        y:          array of shape nTime x nStates
                    observations of the states
        experiment: FGPGM.Experiment object
                    model of the experiment of which the observations have been
                    taken
        nODEParams: scalar
                    number of parameters that need to be fit for the ODE
        gamma:      scalar or vector of length nStates
                    std of the noise matching GP derivatives to those
                    calculated by the ODE
        normalize:  boolean
                    indicates if the kernel hyperparameters have
                    been found using the mean corrected observations
        standardize:    boolean
                    indicates if the kernel hyperparameters have
                    been found using the standardized observations
        
        """
        self.nODEParams = nODEParams
        if np.asarray(gamma).reshape([-1, 1]).size == 1:
            self.gammas = gamma*np.ones(len(kernels))
        else:
            self.gammas = np.asarray(gamma)
        self.CPhis = []
        self.CDashs = []
        self.DashCs = []
        self.CDoubleDashs = []
        self.kernels = kernels
        self.obsNoises = []
        self.time = time
        for kernelID in np.arange(len(kernels)):
            currentKernel = kernels[kernelID]
            self.CPhis.append(currentKernel.getCPhi(time))
            self.CDashs.append(currentKernel.getCPhiDash(time))
            self.DashCs.append(currentKernel.getDashCPhi(time))
            self.CDoubleDashs.append(currentKernel.getCPhiDoubleDash(time))
            self.obsNoises.append(currentKernel.sigma)
        self.As = getAs(self.CDashs, self.DashCs,
                        self.CPhis, self.CDoubleDashs)
        self.Ds = getDs(self.DashCs, self.CPhis)
        self.Lambdas = getLambdaStars(self.gammas, self.time.shape[0])

        for i in np.arange(len(self.As)):
            minEigenval = np.min(np.linalg.eig(self.As[i] + self.Lambdas[i])[0])
            if minEigenval < 0:
                raise Exception(
                    "A + Lambda* is not PSD. Minimum eigenvalue: {}\n".format(
                        minEigenval))


        self.y = y
        self.unfoldedY = np.squeeze(y.reshape([-1, 1], order='F'))
        
        self.experiment = experiment

        self.normalize = normalize
        self.standardize = standardize

        def getDensity(newStates, newParams):
            # new states has already been transformed
            # y is still original
            return calculateLogDensity(self.unfoldedY, newStates,
                                       self.CPhis, self.obsNoises,
                                       self.experiment.f, newParams,
                                       self.As, self.Ds, self.Lambdas,
                                       mean=self.mean,
                                       std=self.std)
        self.getDensity = getDensity
        
        # create mean and std for later use        
        if self.standardize:
            nTime = self.As[0].shape[0]
            nStates = len(self.As)
            assert nStates*nTime == self.y.size
            yMatrix = self.y
            # calculate mean and std
            means = np.mean(yMatrix, axis=0)
            stds = np.std(yMatrix, axis=0)
            assert means.size == nStates
            assert stds.size == nStates
            # reshape mean and std
            self.mean = np.zeros_like(self.unfoldedY)
            self.std = np.zeros_like(self.unfoldedY)
            for i in np.arange(nStates):
                self.mean[nTime*i:nTime*(i+1)] = means[i]*np.ones(nTime)
                self.std[nTime*i:nTime*(i+1)] = stds[i]*np.ones(nTime)
        elif self.normalize:
            nTime = self.As[0].shape[0]
            nStates = len(self.As)
            assert nStates*nTime == self.y.size
            yMatrix = self.y
            # calculate mean
            means = np.mean(yMatrix, axis=0)
            assert means.size == nStates
            self.mean = np.zeros_like(self.unfoldedY)
            for i in np.arange(nStates):
                self.mean[nTime*i:nTime*(i+1)] = means[i]*np.ones(nTime)
            self.std = np.ones_like(self.unfoldedY)
        else:
            self.mean = np.zeros_like(self.unfoldedY)
            self.std = np.ones_like(self.unfoldedY)

    def getFGPGMResults(self, GPPosteriorInit=True, blockNegStates=False,
                        debug=False, theta0=None, thetaMagnitudes=None,
                        nSamples=100000, nBurnin=5000, propStds=None):
        """
        Calculates the optimal ODE parameters and states using the FGPGM
        algorithm
        
        Parameters
        ----------
        GPPosteriorInit:    boolean
                            True: The optimization of the states will be
                            initialized with the GP posterior
                            False: The optimization of the states will be
                            initialized with the observations
        blockNegStates:     boolean
                            if True, bounds on states will be adapted such that
                            even after standardization, they will not become
                            negative
                            Nothing happens if False
        debug:              boolean
                            if True, a plot will be shown comparing the GP
                            posterior for initialization with the observations
                            Nothing happens if False
        theta0:             None or vector of length nODEParams
                            for debugging only.
                            If it is none, the initial guess for theta will
                            be random. If it is not none, the initial guess
                            will be theta0
        thetaMagnitudes:    vector of length nODEParams or None
                            if None, will be set to all ones
                            scaling for the optimization problem. The theta
                            that is fed into the density will be the theta that
                            is fed by the optimizer times 10**thetaMagnitudes
        Returns
        ----------
        newStates:  array of the same shape as self.y
                    inferred states. Direct result of the optimization
        newParams:  vector of length self.nODEParams
                    inferred ODE parameters. Direct result of the optimizaton
        """
        if thetaMagnitudes is None:
            thetaMagnitudes = np.zeros(self.nODEParams)
        thetaMagnitudes = np.asarray(thetaMagnitudes)
        assert thetaMagnitudes.size == self.nODEParams

        def getDensity(optVector):
            """
            calculates the log density, taking just the vector as optimization
            input optVector = [unfoldedStates; parameters], where only the
            states starting from time 1 are flexible
            """
            newStates = optVector[:-self.nODEParams]
            newParams = optVector[-self.nODEParams:]*(10**thetaMagnitudes)
            return -2*newStates.size*self.getDensity(newStates, newParams)

        totalLength = self.y.size + self.nODEParams
        x0 = np.ones(totalLength)

        # initialize theta0
        if theta0 is None:
            theta0 = np.abs(np.random.randn(self.nODEParams))
        else:
            theta0 = theta0 / (10**thetaMagnitudes)
        x0[-self.nODEParams:] = theta0

        np.savetxt("initialGuessForParams.csv", theta0)

        # initialize xInit to observations or with GP posterior
        fullMean = np.zeros_like(self.y)
        # standardize y
        meanMatrix = self.mean.reshape(self.y.shape, order='F')
        stdMatrix = self.std.reshape(self.y.shape, order='F')
        stdY = (self.y - meanMatrix) / stdMatrix
        if GPPosteriorInit:
            print("initialize with GP posterior")
            for i in np.arange(len(self.kernels)):
                currentCPhi = self.kernels[i].getCPhi(self.time)
                currentSigma = self.kernels[i].sigma
                currentY = stdY[:, i]
                GPMean = np.linalg.solve(
                    (currentCPhi + currentSigma*np.eye(currentCPhi.shape[0])),
                    np.dot(currentCPhi, currentY)
                     )
                fullMean[:, i] = GPMean
            # unfold and create as needed by optimizer
            x0[:fullMean.size] = np.squeeze(fullMean.reshape([-1, 1], order='F'))
            np.savetxt("GPPosteriorInit.csv", x0[:fullMean.size])
            np.savetxt("stdMatrix.csv", stdMatrix)
            np.savetxt("meanMatrix.csv", meanMatrix)
            if debug:
                from matplotlib import pyplot as plt
                plt.figure()
                plt.scatter(np.arange(stdY.size),
                            x0[:fullMean.size], marker='x')
                plt.scatter(np.arange(stdY.size),
                         stdY.reshape([-1, 1], order='F'), marker = '.')
                plt.legend(['GPPosterior', 'observations'])
                plt.savefig("./PosteriorInit.png")
                plt.show()
                plt.close()
        else:
            print("initialize states with observations")
            flatObs = np.squeeze(self.y.reshape([-1, 1]), order='F')
            x0[:flatObs.size] = flatObs

        # convert bounds on all states to necessary ones
        xmin, xmax = self.experiment.getBounds(x0.size-self.nODEParams,
                                               self.nODEParams, x0=x0)
        if blockNegStates:
            # adapt bounds such that no negative states will occur
            newZeros = 1e-10*np.ones_like(self.mean)
            newZeros = (newZeros - self.mean) / self.std
            blockCount = 0
            for i in np.arange(x0.size-self.nODEParams):
                if xmin[i] < newZeros[i]:
                    xmin[i] = newZeros[i]
                    blockCount += 1
            if debug:
                print("FGPGM blocked {} potentially negative states".format(blockCount))

        bounds = [(low, high) for low, high in zip(xmin, xmax)]

        # check for legal initialization
        initCorrCount = 0
        for i in np.arange(x0.size):
            if x0[i] < bounds[i][0]:
                x0[i] = bounds[i][0] + bounds[i][1]*1e-5
                initCorrCount += 1
            if x0[i] > bounds[i][1]:
                print("x0[{}] too big".format(i))
        
        print("{} states needed correction, as they were too close to zero".format(initCorrCount))

        if propStds is None:
            print("using standard proposal stds")
            propStds = np.ones(x0.size)
        assert propStds.size == x0.size

        MCMCSamples, nAccepted, nRejected = MCMCWithBounds(
            logTarget=getDensity,
            xInit=x0,
            proposalStds=propStds,
            lowerBounds=xmin,
            upperBounds=xmax,
            nSamples=nSamples,
            nBurnin=nBurnin)
        inferredStuff = np.mean(MCMCSamples, axis=0)

        newStates = inferredStuff[:-self.nODEParams].reshape(
            self.y.shape, order='F')
        # destandardize states
        newStates = stdMatrix*newStates + meanMatrix

        np.savetxt("MCMCMatrix.csv", MCMCSamples)

        newParams = inferredStuff[-self.nODEParams:]*(10**thetaMagnitudes)

        np.savetxt("thetaMagnitudes.csv", thetaMagnitudes)

        stateAccepted = np.asarray(nAccepted[:-self.nODEParams], dtype=np.float)
        stateRejected = nRejected[:-self.nODEParams]
        paramAccepted = np.asarray(nAccepted[-self.nODEParams:], dtype=np.float)
        paramRejected = nRejected[-self.nODEParams:]
        stateRatio = stateAccepted / (stateRejected + stateAccepted)
        paramRatio = paramAccepted / (paramRejected + paramAccepted)
        print("\nstate acceptance mean with std: \n{} +- {}".format(
            np.mean(stateRatio), np.std(stateRatio)))
        print("\nstate acceptance range: \n{} to {}".format(
            np.min(stateRatio), np.max(stateRatio)))
        print("\nparam acceptance mean with std: \n{} +- {}".format(
            np.mean(paramRatio), np.std(paramRatio)))
        print("\nparam acceptance range: \n{} to {}".format(
            np.min(paramRatio), np.max(paramRatio)))

        return newStates, newParams