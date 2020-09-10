import numpy as np
from scipy.spatial import distance_matrix
from pymagi import ArmaVector, ArmaMatrix, ArmaCube, OdeSystem, solveMagiPy
import unittest
from arma import vector, matrix


class SolveMagiTest(unittest.TestCase):
    def test_solve_fn(self):
        fn_system = OdeSystem()
        def fOde(theta, x):
            theta = vector(theta)
            x = matrix(x)

            V = x[:, 0]
            R = x[:, 1]

            Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
            Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
            result = np.stack([Vdt, Rdt], axis=1)

            return ArmaMatrix(result.T.copy())

        def fOdeDx(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
            V = x[:, 0]
            R = x[:, 1]

            resultDx[:,0,0] = theta[2] * (1 - np.square(V))
            resultDx[:,1,0] = theta[2]
            resultDx[:,0,1] = -1.0 / theta[2]
            resultDx[:,1,1] = -1.0*theta[1]/theta[2]
            return ArmaCube(resultDx.T.copy())

        def fOdeDtheta(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
            V = x[:, 0]
            R = x[:, 1]


            resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
            resultDtheta[:,0,1] = 1.0 / theta[2]
            resultDtheta[:,1,1] = -R / theta[2]
            resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)

            return ArmaCube(resultDtheta.T.copy())

        fn_system.fOde = fOde
        fn_system.fOdeDx = fOdeDx
        fn_system.fOdeDtheta = fOdeDtheta
        fn_system.thetaLowerBound = ArmaVector(np.array([0,0,0]))
        fn_system.thetaUpperBound = ArmaVector(np.array([np.inf,np.inf,np.inf]))
        fn_system.name = "FN"
        fn_system.thetaSize = 3


        ydataV = [-0.86, -0.26, 2.14, 1.94, 1.63, 1.75,
                  1.92, 1.39, 1.29, 1.59, 0.63, 0.78, -1.59, -1.92, -1.56, -1.58,
                  -1.26, -1.34, -0.62, -0.39, 1.58, 2.29, 1.69, 1.61, 1.88, 1.57,
                  1.28, 1.09, 1.21, 0.1, -1.66, -2.05, -1.55, -1.81, -1.72, -0.98,
                  -0.77, -0.09, 1.87, 2.18, 1.67]
        ydataR = [0.94, 0.87, 0.62, 0.44,
                  0.07, 0.02, -0.55, -0.09, -0.66, -0.73, -0.73, -0.63, -0.85,
                  -0.55, 0.01, 0.43, 0.4, 0.57, 0.64, 1.26, 1.09, 0.46, 0.13, 0.14,
                  -0.3, -0.53, -0.5, -0.35, -1.03, -1.02, -0.6, -0.61, -0.05, 0.31,
                  0.82, 0.85, 0.64, 1.31, 0.78, 0.47, 0.35]
        ydata = np.stack([np.array(ydataV), np.array(ydataR)], axis=1)
        tvecFull = np.linspace(0, 20, num=81)
        yFull = np.ndarray([81, 2])
        yFull.fill(np.nan)
        yFull[np.linspace(0, 80, num=41).astype(int),:] = ydata
        result_solver = solveMagiPy(
            yFull=ArmaMatrix(yFull).t(),
            odeModel=fn_system,
            tvecFull=ArmaVector(tvecFull),
            sigmaExogenous=ArmaVector(np.ndarray(0)),
            phiExogenous = ArmaMatrix(np.ndarray([0, 0])),
            xInitExogenous = ArmaMatrix(np.ndarray([0, 0])),
            thetaInitExogenous = ArmaVector(np.ndarray(0)),
            muExogenous = ArmaMatrix(np.ndarray([0, 0])),
            dotmuExogenous = ArmaMatrix(np.ndarray([0, 0])),
            priorTemperatureLevel = 1.0,
            priorTemperatureDeriv = 1.0,
            priorTemperatureObs = 1.0,
            kernel = "generalMatern",
            nstepsHmc = 100,
            burninRatioHmc = 0.5,
            niterHmc = 1000,
            stepSizeFactorHmc = 0.1,
            nEpoch = 1,
            bandSize = 20,
            useFrequencyBasedPrior = True,
            useBand = True,
            useMean = False,
            useScalerSigma = False,
            useFixedSigma = False,
            verbose = True)
        phiUsed = matrix(result_solver.phiAllDimensions)
        samplesCpp = matrix(result_solver.llikxthetasigmaSamples.slice(0))


    def test_solve_hes1log(self):
        hes1_system = OdeSystem()
        def fOde(theta, x):
            theta = vector(theta)
            x = matrix(x)

            P = np.exp(x[:, 0])
            M = np.exp(x[:, 1])
            H = np.exp(x[:, 2])

            PMHdt = np.zeros(shape=np.shape(x))
            PMHdt[:,0] = -theta[0]*H + theta[1]*M/P - theta[2]
            PMHdt[:,1] = -theta[3] + theta[4]/(1+np.square(P))/M
            PMHdt[:,2] = -theta[0]*P + theta[5]/(1+np.square(P))/H - theta[6]

            return ArmaMatrix(PMHdt.T.copy())

        def fOdeDx(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
            P = (x[:, 0])
            M = (x[:, 1])
            H = (x[:, 2])

            expMminusP = np.exp(M-P)
            dP = -np.power(1+np.exp(2*P), -2) * np.exp(2*P)*2

            resultDx[:,0,0] = -theta[1]*expMminusP
            resultDx[:,1,0] = theta[1]*expMminusP
            resultDx[:,2,0] = -theta[0]*np.exp(H)
            
            resultDx[:,0,1] = theta[4]*np.multiply(np.exp(-M),dP)
            resultDx[:,1,1] = -theta[4]*np.exp(-M)/(1+np.exp(2*P))
            
            resultDx[:,0,2] = -theta[0]*np.exp(P) + theta[5]*np.multiply(np.exp(-H),dP)
            resultDx[:,2,2] = -theta[5]*np.exp(-H)/(1+np.exp(2*P))

            return ArmaCube(resultDx.T.copy())

        def fOdeDtheta(theta, x):
            theta = vector(theta)
            x = matrix(x)

            resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
            P = (x[:, 0])
            M = (x[:, 1])
            H = (x[:, 2])

            resultDtheta[:,0,0] = -np.exp(H)
            resultDtheta[:,1,0] = np.exp(M-P)
            resultDtheta[:,2,0].fill(-1)

            resultDtheta[:,3,1].fill(-1)
            resultDtheta[:,4,1] = np.exp(-M)/(1+np.exp(2*P))

            resultDtheta[:,0,2] = -np.exp(P)
            resultDtheta[:,5,2] = np.exp(-H)/(1+np.exp(2*P))
            resultDtheta[:,6,2].fill(-1)

            return ArmaCube(resultDtheta.T.copy())

        hes1_system.fOde = fOde
        hes1_system.fOdeDx = fOdeDx
        hes1_system.fOdeDtheta = fOdeDtheta
        hes1_system.thetaLowerBound = ArmaVector(np.array([0,0,0,0,0,0,0]))
        hes1_system.thetaUpperBound = ArmaVector(np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))
        hes1_system.name = "Hes1-log"
        hes1_system.thetaSize = 7


        ydataP = [0.74, 0.78, 1.86, 1.86, 2.2, 1.93, 1.47, 1.03, 0.36,
                  0.88, 1.68, 1.97, 2.15, 1.85, 1.8, 1.47, 0.71]
        ydataM = [0.91, 0.82, 0.71, -0.11, 0.08, -0.45, -0.05, 0.2,
                  0.88, 1.09, 0.3, 0.35, 0.25, -0.23, -0.51, -0.09]
        yFull = np.ndarray([33, 3])
        yFull.fill(np.nan)
        yFull[np.linspace(0, 32, num=17).astype(int),0] = ydataP
        yFull[np.linspace(1, 31, num=16).astype(int),1] = ydataM

        tvecFull = np.linspace(0, 240, num=33)
        result_solver = solveMagiPy(
            yFull=ArmaMatrix(yFull).t(),
            odeModel=hes1_system,
            tvecFull=ArmaVector(tvecFull),
            sigmaExogenous=ArmaVector(np.array([0.15,0.15,0.15])),
            phiExogenous = ArmaMatrix(np.ndarray([0, 0])),
            xInitExogenous = ArmaMatrix(np.ndarray([0, 0])),
            thetaInitExogenous = ArmaVector(np.ndarray(0)),
            muExogenous = ArmaMatrix(np.ndarray([0, 0])),
            dotmuExogenous = ArmaMatrix(np.ndarray([0, 0])),
            priorTemperatureLevel = 1.0,
            priorTemperatureDeriv = 1.0,
            priorTemperatureObs = 1.0,
            kernel = "generalMatern",
            nstepsHmc = 100,
            burninRatioHmc = 0.5,
            niterHmc = 1000,
            stepSizeFactorHmc = 0.1,
            nEpoch = 1,
            bandSize = 20,
            useFrequencyBasedPrior = True,
            useBand = True,
            useMean = False,
            useScalerSigma = False,
            useFixedSigma = True,
            verbose = True)
        phiUsed = matrix(result_solver.phiAllDimensions)
        samplesCpp = matrix(result_solver.llikxthetasigmaSamples.slice(0))

