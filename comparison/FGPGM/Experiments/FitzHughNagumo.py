# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
implements the FitzHugh Nagumo dynamics
"""
from ..Experiment import Experiment
import numpy as np

class FHN(Experiment):
    def __init__(self):
        pass
    
    def f(self, x, theta):
        """
        \dot{V} = \phi*(V - V**3/3 - R)
        \dot{R} = 1/\phi * (V- \alpha + \beta*R)
        will return derivatives in a vector
        """
        V = x[0]
        R = x[1]
        alpha = theta[0]
        beta = theta[1]
        phi = theta[2]
        firstODE = phi*(V - V**3/3 + R)
        secondODE = -1/phi*(V - alpha + beta*R)
        return np.asarray([firstODE, secondODE])

    def getBounds(self, nStates, nParams, x0=None):
        # mby move to main optimization body??
        """
        returns a list of bounds for states and parameters to constrain
        optimization
        """
        xmin = []
        xmax = []
        for i in np.arange(nStates):
            xmin.append(-3)
            xmax.append(3)
        for i in np.arange(nParams):
            xmin.append(1e-4)
            xmax.append(10)
        return xmin, xmax