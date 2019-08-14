import numpy as np
from scipy.spatial import distance_matrix
from pygpds import ArmaVector, ArmaMatrix, OdeSystem
import unittest
from arma import vector, matrix


class OdeSystemTest(unittest.TestCase):
    def test_constructor(self):
        ode_system = OdeSystem()
        def fOde(theta, x):
            theta = vector(theta)
            x = matrix(x)
            return theta * x
        ode_system.fOde = fOde
