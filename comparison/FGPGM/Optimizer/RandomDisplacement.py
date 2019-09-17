# -*- coding: utf-8 -*-
# Author: Philippe Wenk <philippewenk@hotmail.com>
# based on scipy _basinhopping.py
import numpy as np

class RandomDisplacement(object):
    """
    Add a random displacement of maximum size, stepsize, to the coordinates
    update x inplace

    Parameters
    ----------
    stepsize:   float, optional
                stepsize to be taken. Will be adapted by the optimizer later.
    bounds:     list of tuples
                bounds on the variables.
                If None, global optimization is performed.
                If [(a,b), (c,d), ...], all bounds need to be specified. Length
                of bounds must be the same as the length of the parameter
                vector.
    """
    def __init__(self, stepsize=100, bounds=None):
        self.stepsize = stepsize
        self.random_state = np.random
        self.bounds = np.asarray(bounds)
        if self.bounds.shape[0] == self.bounds.size:
            if not self.bounds.shape[0]==2:
                raise TypeError("incorrect data type chosen for bounds.")
        else:
            if not self.bounds.shape[1]==2:
                raise TypeError("incorrect data type chosen for bounds.")
        self.min = self.bounds[:, 0]
        self.max = self.bounds[:, 1]


    def __call__(self, x):
        # check if x one dimensional
        if not (self.bounds.shape[0] == x.size or \
            (x.size==1 and self.bounds.size==2)):
            raise ValueError(
                "incorrect shape of data or bounds, {} != {}".format(
                    self.bounds.shape[0], x.size))
        newX = np.zeros_like(x)
        # iterate through all components. check each
        for i in np.arange(x.size):
            currentStepSize = self.stepsize
            counter = 0
            success = False
            #decrease step size until valid step found. Do this 10 times
            #if all fails, use random number uniformly distributed within bound
            while(counter < 10 and not success):
                for j in np.arange(10):
                    newX[i] = x[i] + self.random_state.uniform(
                        -currentStepSize, currentStepSize, 1)
                    # exit if new value within steps
                    if newX[i] < self.max[i] and newX[i] > self.min[i]:
                        success = True
                        break
                counter += 1
                currentStepSize = float(currentStepSize)/2
            if not success:
                newX[i] = self.random_state.uniform(
                    self.min[i], self.max[i], 1)
                print("Hacked")
        return newX