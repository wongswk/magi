# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
"""
Implements the Protein Transduction dynamics
"""
from ..Experiment import Experiment
import numpy as np

class ProteinTransduction(Experiment):
    def __init__(self):
        pass
    
    def f(self, x, theta):
        """
        x: (S, S_d, R, RS, Rpp)
        theta: (k1, k2, k3, k4, V, Km)
        
        Creates the odes for given states and thetas
        """
        # make parameters human readable
        k1 = theta[0]
        k2 = theta[1]
        k3 = theta[2]
        k4 = theta[3]
        V = theta[4]
        Km = theta[5]
        # make states human readable
        S = x[0]
        R = x[2]
        RS = x[3]
        Rpp = x[4]
        # create ODEs
        dS = -k1*S - k2*S*R + k3*RS
        dS_d = k1*S
        dR = -k2*S*R + k3*RS + float(V*Rpp)/(Km + Rpp + 1e-6)
        dRS = k2*S*R - k3*RS - k4*RS
        dRpp = k4*RS - float(V*Rpp)/(Km + Rpp + 1e-6)
        return np.asarray([dS, dS_d, dR, dRS, dRpp])

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
            xmin.append(-5)
            xmax.append(5)
        for i in np.arange(nParams):
            xmin.append(0)
            xmax.append(4)

        return xmin, xmax
