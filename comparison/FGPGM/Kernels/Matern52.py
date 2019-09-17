# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Implements the necessary functions of the Kernel class for the RBF kernel.
"""

import numpy as np
from ..Kernel import Kernel

class Matern52(Kernel):
    def __init__(self, theta=None, sigma=None, nugget=1e-4):
        """
        Parameters
        ----------
        theta:  vector of shape n_hyperparameters
                vector containing the hyperparameters of this kernel
                theta[0] is the multiplicative constant sigma_f**2
                theta[1] is the lengthscale l
                Both will be guaranteed to be positive during optimization
        sigma:  scalar
                std of the observation noise
        nugget: scalar
                small float that is added to the diagonal of the kernel matrix
                CPhi to guarantee positive definiteness also numerically.
        """
        self.theta = theta
        self.sigma = sigma
        self.nugget = nugget

    def getNTheta(self):
        """
        Returns the amount of parameters needed by this kernel in the theta
        vector
        """
        return 2

    def k(self, time1, time2):
        """
        returns the correlation between time1 and time2 for the specific kernel
        this does not yet add the observation noise
        """
        time1 = float(time1)
        time2 = float(time2)
        r = np.abs(time1 - time2)
        return self.theta[0]\
             * (1 + np.sqrt(5)/self.theta[1]*r + 5./3./self.theta[1]**2*r**2) \
             * np.exp(-np.sqrt(5)/self.theta[1]*r)
        
    def CDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time2, used in the C_Phi' matrix
        """
        return  self.DashC(time2, time1)
    
    def DashC(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time1, used in the 'C_Phi matrix
        """
        time1 = float(time1)
        time2 = float(time2)
        r = np.abs(time1 - time2)
        return - self.theta[0]*5./3./self.theta[1]**2*(time1 - time2) \
               * (1 + np.sqrt(5)/self.theta[1]*r) \
               * np.exp(-np.sqrt(5)/self.theta[1]*r)

    def CDoubleDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to both times, used in the C_Phi'' matrix
        """
        time1 = float(time1)
        time2 = float(time2)
        r = np.abs(time1 - time2)
        return   self.theta[0]*5./3./self.theta[1]**2 \
               * (1 + np.sqrt(5)/self.theta[1]*r - 5./self.theta[1]**2*r**2) \
               * np.exp(-np.sqrt(5)/self.theta[1]*r)

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
        # TODO: justify bounds theoretically
        upperBoundSigmaF = (np.max(y) - np.min(y))**2
        upperBoundLengthscale = time[1]*100 #TODO: do 10 again?
        upperBoundStd = np.max(y) - np.min(y)
        lowerBoundLengthscale = time[1]
        bounds = [(1e-4, upperBoundSigmaF),
                  (lowerBoundLengthscale, upperBoundLengthscale),
                  (1e-3, upperBoundStd)
                  ]
        return bounds
