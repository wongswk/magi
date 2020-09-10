import numpy as np
from scipy.spatial import distance_matrix
from pymagi import ArmaVector, ArmaMatrix, basic_hmcC, lp, hmcstate
import unittest
from arma import vector, matrix


class HmcTest(unittest.TestCase):
    def test_hmc(self):
        def llk(x):
            x = vector(x)
            llik_value = -0.5 * np.sum(np.square(x))
            llik_gradient = -x
            llik = lp()
            llik.value = llik_value
            llik.gradient = ArmaVector(llik_gradient)
            return llik
        out = basic_hmcC(lpr=llk,
                         initial=ArmaVector([0.2, 0.2]),
                         step=ArmaVector([0.1, 0.1]),
                         lb=ArmaVector([-10, -10]),
                         ub=ArmaVector([np.Inf, np.Inf]),
                         nsteps=200,
                         traj=False)
        self.assertEqual(out.final.size(), 2)
        self.assertEqual(llk(out.final).value, out.lprvalue)
