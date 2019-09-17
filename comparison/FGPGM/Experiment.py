# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>

from scipy.integrate import odeint
import numpy as np
import os
from matplotlib import pyplot as plt

class Experiment(object):
    def __init__(self):
        pass
    
    def f(self, x, theta):
        """
        representing the ODEs. Will return a vector with same shape as the
        state vector x, containing the derivatives of each state w.r.t. time
        for the given ODE parameters theta
        """
        raise NotImplementedError(
            "f has not been implemented for this Experiment")

    def getConstraints(self, nStates, nParams):
        # currently untested and unused
        """
        returns a dict of functions that receive a vector with flattened
        states and parameters, representing constraints as required by the
        optimizer
        (See: https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#constrained-minimization-of-multivariate-scalar-functions-minimize)
        
        Input vector of the functions will be a vector [unfoldedX, theta]
        where theta are the ODE parameters as taken by theta and unfoldedX are
        the states unfolded like this:
        unfoldedX = [x1[t0], x1[t1], ..., x1[tEnd], x2[t0]...]
        """
        raise NotImplementedError(
            "constraints have not been implemented for this Experiment")

    def getBounds(self, nStates, nParams):
        """
        returns a list of bounds for states and parameters to constrain
        optimization
        """
        raise NotImplementedError(
            "bounds have not been implemented for this Experiment")

#    def sampleTrajectory(self, XInit, tEnd, dt, theta, obsNoiseStd):
#        """
#        Creates a trajectory using a numerical integrator
#        
#        Parameters
#        ----------
#        XInit:          vector of length nStates
#                        initial values of the states at time zero
#        tEnd:           scalar
#                        end time of the experiment
#        dt:             scalar
#                        time steps at which an observation should be sampled
#        theta:          vector of length nParameters
#                        parameters for the ODE
#        obsNoiseStd:    scalar
#                        std of the noise on observations
#        Returns
#        ----------
#        x:  array of shape nTime x nStates
#            true states as returned by the integrator
#        y:  array of shape nTime x nStates
#            noisy observations of the true states
#        """
#        def fODE(x, time):
#            return self.f(x, theta)
#        time = np.arange(0, tEnd+0.5*dt, dt)
#        x = odeint(fODE, XInit, time, rtol=1e-8, mxstep=5000000) # huge for stiff problems
#        noise = np.random.randn(x.shape[0], x.shape[1])
#        noise = noise*obsNoiseStd
#        y = x + noise
#        return x, y
    
    def sampleTrajectoryNonUniform(self, XInit, theta, time, obsNoiseStd=None,
                                   SNR=None, plotting=None):
        """
        Creates a trajectory using a numerical integrator
        
        Pararmeters
        -----------
        XInit:          vector of length nStates
                        initial values of the states at time zero
        theta:          scalar
                        end time of the experiment
        time:           vector
                        time points at which the trajectory should be observed
        obsNoiseStd:    scalar or None
                        std of the noise on observations
                        if None, SNR will be used instead
        SNR:            scalar or None
                        if obsNoisStd is not None, this will be ignored
                        Else, this specifies the signal to noise ratio for the
                        observations y
        plotting:       None or string
                        path at which the experiment plots should be stored
                        if None, no plots will be created.
        Returns
        ----------
        x:  array of shape nTime x nStates
            true states as returned by the integrator
        y:  array of shape nTime x nStates
            noisy observations of the true states
        """
        def fODE(x, time):
            return self.f(x, theta)
        x = odeint(fODE, XInit, time, rtol=1e-8, mxstep=5000000) # huge for stiff problems
        noise = np.random.randn(x.shape[0], x.shape[1])
        if obsNoiseStd is None:
            signalStds = np.std(x, axis=0)
            obsNoiseStds = signalStds/np.sqrt(SNR)
            obsNoiseStds = obsNoiseStds.reshape([1, -1])
            obsNoiseStds = np.repeat(obsNoiseStds, x.shape[0], axis=0)
        else:
            obsNoiseStds = np.ones_like(x)*obsNoiseStd
        noise = noise*obsNoiseStds
        y = x + noise
        if plotting is not None:
            if not os.path.exists(plotting):
                os.makedirs(plotting)
            for state in np.arange(x.shape[1]):
                # set ticks
                fig = plt.figure()
                plot = fig.add_subplot(111)
                # set cross thickness
                plot.scatter(time, y[:, state], c='k', marker='.', s=100,
                            linewidths=1)
                # set line thickness
                plot.plot(time, x[:, state], 'r', linewidth=2)
                # set label fontsize
                plt.xlabel("time",
                           fontsize=20)
                plt.ylabel("state {}".format(state+1),
                           fontsize=20)
                plt.setp(plot.get_xticklabels(),
                         fontsize=20)
                plt.setp(plot.get_yticklabels(),
                         fontsize=20)
                plt.tight_layout()
                plt.savefig(os.path.join(plotting,
                                         "state{}.png".format(state)),
                            dpi=300)
                
                plt.close()
        return x, y