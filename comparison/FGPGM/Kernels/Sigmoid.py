# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Implements the necessary functions of the Kernel class for the Sigmoid kernel.
"""
import numpy as np
from ..Kernel import Kernel

class Sigmoid(Kernel):
    def __init__(self, theta=None, sigma=None, nugget=1e-4):
        """
        Sigmoid kernel calculating the covariance function
        
        k(t1, t2) = sigmaF**2 * arcsin((a+b*t1*t2) / 
                    sqrt[(a + b*t1**2 + 1)(a + b*t2*t2 + 1)]

        Parameters
        ----------
        theta:  vector of shape n_hyperparameters
                vector containing the hyperparameters of this kernel
                theta[0] is sigmaF**2
                theta[1] is a
                theta[2] is b
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
        return 3

    def k(self, time1, time2):
        """
        returns the correlation between time1 and time2 for the specific kernel
        this does not yet add the observation noise
        """
        time1 = float(time1)
        time2 = float(time2)
        # write parameters in humand readable form
        sigmaFSq = self.theta[0]
        return sigmaFSq*np.arcsin(self._getZ(time1, time2))

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
        return [(1e-4, 100),
                 (1e-4, 250),
                 (1e-4, 250),
                 (1e-4, 2)]

    def CDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time2, used in the C_Phi' matrix
        """
        sigmaFSq = self.theta[0]
        return sigmaFSq / \
               np.sqrt(1-self._getZ(time1, time2)**2) * \
               self._getDZDt2(time1, time2)
    
    def DashC(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to time1, used in the 'C_Phi matrix
        """
        sigmaFSq = self.theta[0]
        return sigmaFSq / \
               np.sqrt(1-self._getZ(time1, time2)**2) * \
               self._getDZDt1(time1, time2)

    def CDoubleDash(self, time1, time2):
        """
        returns the derivative of the correlation between time1 and time2 with
        respect to both times, used in the C_Phi'' matrix
        """
        sigmaFSq = self.theta[0]
        return sigmaFSq / \
               np.sqrt(1-self._getZ(time1, time2)**2) * \
               (self._getZ(time1, time2) / \
                (1-self._getZ(time1, time2)**2) * \
                self._getDZDt1(time1, time2)*self._getDZDt2(time1, time2)
                +
                self._getDZDt1Dt2(time1, time2)
               )
    
    def _getZ(self, time1, time2):
        """
        calculates the argument of the arcsin
        """
        a = self.theta[1]
        b = self.theta[2]
        return (a + b*time1*time2)/self._getZNorm(time1, time2)

    def _getZNorm(self, time1, time2):
        """
        calculates the denominator of the argument of the arcsin
        """
        a = self.theta[1]
        b = self.theta[2]
        return np.sqrt((a + b*time1**2 + 1)*(a + b*time2**2 + 1))

    def _getDZDt1(self, time1, time2):
        """derivative of Z w.r.t. time1"""
        a = self.theta[1]
        b = self.theta[2]
        firstSummand = b*time2/self._getZNorm(time1, time2)
        secondSummand = b*time1*self._getZ(time1, time2)/(a + b*time1**2 + 1)
        return firstSummand - secondSummand

    def _getDZDt2(self, time1, time2):
        """derivative of Z w.r.t. time2"""
        a = self.theta[1]
        b = self.theta[2]
        firstSummand = b*time1/self._getZNorm(time1, time2)
        secondSummand = b*time2*self._getZ(time1, time2)/(a + b*time2**2 + 1)
        return firstSummand - secondSummand

    def _getDZDt1Dt2(self, time1, time2):
        """second order derivative of Z w.r.t. time1 and time2"""
        a = self.theta[1]
        b = self.theta[2]
        firstSummand = b / self._getZNorm(time1, time2)
        secondSummand = -b**2*time2**2/(
            self._getZNorm(time1, time2)*(a + b*time2**2 + 1))
        thirdSummand = -b*time1/(a + b*time1**2 + 1) * \
                       self._getDZDt2(time1, time2)
        return firstSummand + secondSummand + thirdSummand
