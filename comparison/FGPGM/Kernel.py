# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Defines the Kernel class. Instances should be implemented in the Kernel folder.
"""

import numpy as np
from scipy.optimize import basinhopping
from .Optimizer.RandomDisplacement import RandomDisplacement

from scipy.misc import derivative

class Kernel(object):
    def __init__(self, theta=None, sigma=None, nugget=1e-4):
        """
        Parameters
        ----------
        theta:  vector of shape n_hyperparameters
                vector containing the hyperparameters of this kernel
        sigma:  scalar
                std of the observation noise
        nugget: scalar
                small float that is added to the diagonal of the kernel matrix
                CPhi to guarantee positive definiteness also numerically.
        """
        self.theta = theta
        self.sigma = sigma
        self.nugget = nugget

    def getNTheta():
        """
        Returns the amount of parameters needed by this kernel in the theta
        vector
        """
        raise NotImplementedError(
            "getNTheta has not been implemented for this kernel.")

    def setHyperparams(self, theta, sigma, newNugget=False):
        """
        Sets the hyperparameters of this kernel
        """
        self.theta = theta
        self.sigma = sigma
        if not (isinstance(newNugget, bool) and not newNugget):
            self.nugget = newNugget

    def negLogLikelihood(self, params, time, normalize, standardize, y):
        """
        calculates the negative Log likelihood needed to maximize the evidence
        of the data for the hyperparameters
        """
        # get C matrix without disturbing current hyperparameters
        sigma = params[-1]
        theta = params[:-1]
        thetaOld = self.theta
        sigmaOld = self.sigma
        self.setHyperparams(theta, sigma)
        C = self.getCPhi(time)
        C = C + sigma**2*np.eye(C.shape[0])
        self.setHyperparams(thetaOld, sigmaOld)
        # calculate actual log likelihood
        if standardize:
            yNorm = (y - np.mean(y)) / np.std(y)
        elif normalize:
            yNorm = (y - np.mean(y))
        else:
            yNorm = y
        sum1 = np.prod(np.linalg.slogdet(C))
        sum2 = np.dot(yNorm, np.linalg.solve(C, yNorm))
        assert sum1 is not np.nan
        assert sum2 is not np.nan
        return (sum1 + sum2) / y.size                     

    def getBounds(self, y, time):
        """
        creates the bounds for the optimization of the hyperparameters.
        
        Parameters
        ----------
        y:          vector
                    observation of the states. Target of the regression
        time:       vector
                    time points of the observations. Input of the regression
        Returns
        ----------
        bounds: list of theta.size + 1 pairs of the form 
                (lowerBound, upperBound), representing the bounds on the 
                kernel hyperparameters in theta, while the last one is the
                bound on sigma
        """
        raise NotImplementedError(
            "getBounds has not yet been implemented for this kernel!"
            )

    def learnHyperparams(self, theta0, sigma0, y, time, normalize=False,
                         standardize=False, T=1, newNugget=False, anneal=False,
                         annealArgs={}, basinIter=100):
        """
        Learns the hyperparameters by maximizing the marginal likelihood of the
        data y

        Parameters
        ----------
        theta0:     vector
                    initial guess for parameters for optimization
        sigma0:     scalar
                    initial guess for noise for optimization
        y:          vector of length nObs or array of shape nObs x nReps
                    observation of the states. Target of the regression
                    if y is an array, it is assumed that the observations
                    come from different, independent experiments on the same
                    time scale.
        time:       vector of length nObs
                    time points of the observations. Input of the regression
        normalize:  boolean
                    if True, hyperparameters will be optimized for the
                    mean corrected observation.
                    if False, hyperparameters will be optimized directly
        standardize:    boolean
                        if True, hyperparameters will be optimized for the
                        standardized observations. normalize will be ignored
                        if False, hyperparameters will be optimized as
                        specified by normalize keyword.
        T:          scalar
                    Temperature for the basinhopping optimization
        plotParams: vector of length 5
                    if not None, this code will create a heatmap, plotting
                    the two parameters against each other with the
                    negLogLikelihood as value.
                    if None, nothing will happen.
                    vector is organized as [xmin, xmax, ymin, ymax, trueNoise]
        newNugget:  False or float
                    if false, the old nugget will be used
                    if float, the old nugget will be overwritten
                    nugget is the small number that is added to the GP prior
                    covariance matrix to guarantee positive definiteness
                    also numerically
        anneal:     bool
        annealArgs: dict
        basinIter:  scalar
                    if no annealing is performed, basinIter iterations of
                    basinhopping will be done instead
        """
        # define optimization target
        if y.size == y.shape[0]:
            def negLogLikelihood(params):
                return self.negLogLikelihood(params, time, normalize,
                                             standardize, y)
        else:
            def negLogLikelihood(params):
                # for multiple trajectories, just add likelihood of each
                # run. Assumes one GP per trajectory and mean likelihood as
                # optimization target
                likelihoods = 0
                for i in np.arange(y.shape[1]):
                    likelihoods += self.negLogLikelihood(
                        params, time, normalize, standardize, y[:, i])
                return likelihoods
        # set optimizer settings
        bounds = self.getBounds(y, time)
                  
        # set nugget
        if not (isinstance(newNugget, bool) and not newNugget):
            self.nugget = newNugget
        else:
            print(newNugget)
        
        print("using L-BFGS-B as hyperparameter optimizer")
        # include method and bounds
        args = dict(method="L-BFGS-B", bounds=bounds)
        # default options from scipy
        options={'disp': None,
                 'maxls': 20,
                 'iprint': -1,
                 'gtol': 1e-05,
                 'eps': 1e-08,
                 'maxiter': 15000,
                 'ftol': 2.220446049250313e-09,
                 'maxcor': 10,
                 'maxfun': 15000}
        # ftol: relative difference in function value accepted for convergence
        options['ftol'] = 2.220446049250313e-09
        # maximum number of function evaluations
        options['maxfun'] = 200000
        # flag to control showing of convergence messages
        args['options'] = options
        x0 = np.zeros(theta0.size + 1)
        x0[:-1] = theta0
        x0[-1:] = sigma0
        if sigma0 < 1e-3:
            sigma0 = 1e-3
        def printAcceptance(x, f, accept):
            if accept:
                print("YES: {} @ {}".format(f, x))
            else:
                print("Nope: {} @ {}".format(f, x))
        
        if not anneal:            
            minimum = basinhopping(negLogLikelihood, x0, T=T,
                                   minimizer_kwargs=args,
                                   take_step=RandomDisplacement(bounds=bounds),
                                   niter=int(basinIter),
                                   callback=printAcceptance
                                   )
        else:
            from .Optimizer.Annealer import simulatedAnneal
            basinArgs = dict()
            basinArgs['minimizer_kwargs'] = args
            basinArgs['take_step'] = RandomDisplacement(bounds=bounds)
            basinArgs['callback'] = "./logData/HyperparameterOptimization/"
            minimums = simulatedAnneal(
                negLogLikelihood, x0,
                Temps=10**np.linspace(2, -5, 8),
                iterations=50,
                basinArgs=basinArgs
                )
            minimum = minimums[-1]
            for i in np.arange(len(minimums)):
                print(minimums[i].x)
        print("Kernel optimization output: ")
        print(minimum)
        print("\n")
        optVector = minimum.x
        self.theta = optVector[:-1]
        self.sigma = optVector[-1]
        # check for positive semidefinite
        C = self.getCPhi(time)
        minEigenvalue = np.min(np.linalg.eig(C)[0])
        print("minimum eigenvalue = {}".format(minEigenvalue))
        if (minEigenvalue < 1e-5):
            print("\n\nRECOMMENDATION: USE BIGGER NUGGET\n\n")
        C = C + (self.sigma**2)*np.eye(C.shape[0])
        try:
            # test for psd
            np.linalg.cholesky(C)
        except:
            print("matrix not positive semidefinite")

    def k(self, time1, time2):
        """
        returns the correlation between time1 and time2 for the specific kernel
        this does not yet add the observation noise
        """
        raise NotImplementedError("k has not been implemented for this kernel")
    
    def CDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time2, used in the C_Phi' matrix
        """
        raise NotImplementedError(
            "CDash has not been implemented for this kernel")
    
    def DashC(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time1, used in the 'C_Phi matrix
        """
        raise NotImplementedError(
            "DashC has not been implemented for this kernel")
    
    def CDoubleDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to both times, used in the C_Phi'' matrix
        """
        raise NotImplementedError(
            "CDoubleDash has not been implemented for this kernel")
    
    def getCPhi(self, time):
        """
        returns the correlation matrix of the GP using this kernel
        """
        C_Phi = np.zeros([time.shape[0], time.shape[0]])
        for i in np.arange(time.shape[0]):
            for j in np.arange(time.shape[0]):
                C_Phi[i, j] = self.k(time[i], time[j])
        return C_Phi + self.nugget*np.eye(time.shape[0])

    def getCPhiDash(self, time):
        """
        returns the derivative of C_Phi w.r.t. the second time argument
        """
        C_PhiDash = np.zeros([time.shape[0], time.shape[0]])
        for i in np.arange(time.shape[0]):
            for j in np.arange(time.shape[0]):
                C_PhiDash[i, j] = self.CDash(time[i], time[j])
        return C_PhiDash

    def getDashCPhi(self, time):
        """
        returns the derivative of C_Phi w.r.t. the first time argument
        """
        DashC_Phi = np.zeros([time.shape[0], time.shape[0]])
        for i in np.arange(time.shape[0]):
            for j in np.arange(time.shape[0]):
                DashC_Phi[i, j] = self.DashC(time[i], time[j])
        return DashC_Phi

    def getCPhiDoubleDash(self, time):
        """
        returns the derivative of C_Phi w.r.t. both time arguments
        """
        C_PhiDoubleDash = np.zeros([time.shape[0], time.shape[0]])
        for i in np.arange(time.shape[0]):
            for j in np.arange(time.shape[0]):
                C_PhiDoubleDash[i, j] = self.CDoubleDash(time[i], time[j])
        return C_PhiDoubleDash
    
    def testValidity(self, verbose=False, dt=1e-2, tol=1e-5):
        """
        function to check for the most common mistakes in creating a kernel.
        Will test kernel parameter existence, symmetry of kernel function
        and compare derivatives with numerical approximations.
        
        Parameters
        ----------
        verbose:    boolean
                    if True, this function will tell you what it is doing.
        dt:         scalar
                    dx argument for scipy.misc.derivative
        tol:        scalar
                    relative tolerance between numerical and analytical gradient
        """

        "Test kernel parameters"
        if self.theta is None:
            raise TypeError(
                "theta is currently None. Please initialize with appropriate values")
        if verbose:
            print("current theta: {}".format(self.theta))
        if self.sigma is None:
            raise TypeError(
                "sigma is currently None. Please initialize with appropriate value.")
        if verbose:
            print("current sigma: {}".format(self.sigma))
            print("current nugget: {}".format(self.nugget))

        "check for symmetry"
        if verbose:
            print("\nStart symmetry testing")
        for i in np.arange(10):
            times = np.abs(3*(np.random.randn(2)))
            if not np.allclose(self.k(times[0], times[1]), 
                               self.k(times[1], times[0])):
                raise StandardError(
                    "kernel function is not symmetfic for times {}, {}".format(
                        times[0], times[1]))
        if verbose:
            print("Successful termination of symmetry testing")

        "Numerical derivative check"
        if verbose:
            print("\nStart derivative testing")
            print("spacing: {}".format(dt))

        def getNumCDash(time1, time2):
            def func(t2):
                return self.k(time1, t2)
            return derivative(func, time2, dx=dt)

        def getNumDashC(time1, time2):
            def func(t1):
                return self.k(t1, time2)
            return derivative(func, time1, dx=dt)

        def getNumCDD(time1, time2):
            def func(t2):
                return getNumDashC(time1, t2)
            return derivative(func, time2, dx=dt)

        for i in np.arange(10):
            times = np.abs(3*np.random.randn(2))
            # test derivative w.r.t. second time argument
            kernelCDash = self.CDash(times[0], times[1])
            numCDash = getNumCDash(times[0], times[1])
            if not np.allclose(kernelCDash, numCDash, rtol=tol):
                raise StandardError(
                    "implemented CDash does not agree with numerical test for" + 
                    " times\n{}, {}\nwith values\n{} != {}".format(
                        times[0], times[1], kernelCDash, numCDash)
                    )
            # test derivative w.r.t. first time argument
            kernelDashC = self.DashC(times[0], times[1])
            numDashC = getNumDashC(times[0], times[1])
            if not np.allclose(kernelDashC, numDashC, rtol=tol):
                raise StandardError(
                    "implemented DashC does not agree with numerical test for" +
                    " times\n{}, {}\nwith values\n{} != {}".format(
                        times[0], times[1], kernelDashC, numDashC)
                    )
            # test derivative w.r.t. both time arguments
            kernelCDD = self.CDoubleDash(times[0], times[1])
            numCDD = getNumCDD(times[0], times[1])
            if not np.allclose(kernelCDD, numCDD, rtol=tol):
                raise StandardError(
                    "implemented CDoubleDash does not agree with numerical" +
                    " test for times\n{}, {}\n with values {} != {}".format(
                        times[0], times[1], kernelCDD, numCDD)
                    )
        if verbose:
            print("Successful temination of derivative testing.")