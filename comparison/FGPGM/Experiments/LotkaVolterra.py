# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
implements the Lotka Volterra dynamics
"""
from ..Experiment import Experiment
import numpy as np

class LotkaVolterra(Experiment):
    def __init__(self):
        pass
    
    def f(self, x, theta):
        """
        \dot{x_0} = x_0 * (theta_0 - theta_1*x_1)
        \dot{x_1} = -x_1 * (theta_2 - theta_3*x_0)
        will return derivatives in a vector
        """
        firstODE = x[0]*(theta[0] - theta[1]*x[1])
        secondODE = -x[1]*(theta[2] - theta[3]*x[0])
        return np.asarray([firstODE, secondODE])

    def getBounds(self, nStates, nParams, x0=None):
        """
        returns a list of bounds for states and parameters to constrain
        optimization. If standardization is used later on, the bounds on the
        states are for the standardized states, i.e. [-3, 3] on a state
        variable with mean 2 and std 3 will mean bounds of [-7, 11] in original
        coordinates.
        """
        xmin = []
        xmax = []
        for i in np.arange(nStates):
            xmin.append(-3)
            xmax.append(3)
        for i in np.arange(nParams):
            xmin.append(0)
            xmax.append(100)
        return xmin, xmax
