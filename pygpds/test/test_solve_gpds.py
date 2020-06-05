import numpy as np
from scipy.spatial import distance_matrix
from pygpds import ArmaVector, ArmaMatrix, ArmaCube, OdeSystem
import unittest
from arma import vector, matrix


class SolveGpdsTest(unittest.TestCase):
    def test_solve_gpds_py(self):
        fn_system = OdeSystem()
        def fOde(theta, x):
            theta = vector(theta)
            x = matrix(x)

            V = x.loc[:, 0]
            R = x.loc[:, 1]

            Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
            Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
            result = np.stack([Vdt, Rdt], axis=1)

            return ArmaMatrix(result)

        def fOdeDx(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
            V = x.loc[:, 0]
            R = x.loc[:, 1]

            resultDx[:,0,0] = theta[2] * (1 - np.square(V))
            resultDx[:,1,0] = theta[2]
            resultDx[:,0,1] = -1.0 / theta[2]
            resultDx[:,1,1] = -1.0*theta[1]/theta[2]
            return ArmaCube(resultDx)

        def fOdeDtheta(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
            V = x.loc[:, 0]
            R = x.loc[:, 1]


            resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
            resultDtheta[:,0,1] = 1.0 / theta[2]
            resultDtheta[:,1,1] = -R / theta[2]
            resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)

            return ArmaCube(resultDtheta)

        fn_system.fOde = fOde
        fn_system.fOdeDx = fOdeDx
        fn_system.fOdeDtheta = fOdeDtheta
        fn_system.thetaLowerBound = ArmaVector(np.array([0,0,0]))
        fn_system.thetaUpperBound = ArmaVector(np.array([np.inf,np.inf,np.inf]))
        fn_system.name = "FN"
        fn_system.thetaSize = 3
