# Replication script for Python

import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# Add path to MAGI library directory (containing pymagi.so) here if necessary
# sys.path.append('...')

from arma import ode_system, testDynamicalModel, setDiscretization, gpsmoothing
from magi import MagiSolver

## Hes1 example

def hes1modelOde(theta, x, tvec):
    P = x[:, 0]
    M = x[:, 1]
    H = x[:, 2]
    
    PMHdt = np.zeros(shape=np.shape(x))
    PMHdt[:,0] = -theta[0]*P*H + theta[1]*M - theta[2]*P
    PMHdt[:,1] = -theta[3]*M + theta[4]/(1+np.square(P))
    PMHdt[:,2] = -theta[0]*P*H + theta[5]/(1+np.square(P)) - theta[6]*H
    return PMHdt
    
true_theta = [0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3]
true_x0 = [1.439, 2.037, 17.904]
true_sigma=[0.15, 0.15, np.NaN]

tvecObs = np.linspace(0, 240, 33)

tvecOde = np.linspace(0, 240, 24001)
sol = solve_ivp(lambda t, y: hes1modelOde(true_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecOde[-1]], y0=true_x0, t_eval=tvecOde, vectorized=True)

np.random.seed(12321)

y = sol.y.transpose()
y = y[np.linspace(0,24000,33).astype(int),:]

y[:, 0] *=  np.exp( np.random.normal(0, true_sigma[0], 33))
y[:, 1] *=  np.exp( np.random.normal(0, true_sigma[1], 33))

y[np.linspace(1, 31, num=16).astype(int), 0] = np.nan
y[np.linspace(0, 32, num=17).astype(int), 1] = np.nan
y[:, 2] = np.nan

plt.plot(tvecOde, sol.y[0, :], label = 'P')
plt.plot(tvecOde, sol.y[1, :], label = 'M')
plt.plot(tvecOde, sol.y[2, :], label = 'H')
plt.xlabel('Time')
plt.ylabel('Level')
plt.legend()
plt.scatter(tvecObs , y[:, 0])
plt.scatter(tvecObs , y[:, 1])
plt.show()

y = np.log(y)

def hes1logmodelOde(theta, x, tvec):
    P = np.exp(x[:, 0])
    M = np.exp(x[:, 1])
    H = np.exp(x[:, 2])
    
    PMHdt = np.zeros(shape=np.shape(x))
    PMHdt[:,0] = -theta[0]*H + theta[1]*M/P - theta[2]
    PMHdt[:,1] = -theta[3] + theta[4]/(1+np.square(P))/M
    PMHdt[:,2] = -theta[0]*P + theta[5]/(1+np.square(P))/H - theta[6]
    return PMHdt


def hes1logmodelDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    logP = (x[:, 0])
    logM = (x[:, 1])
    logH = (x[:, 2])

    expMminusP = np.exp(logM-logP)
    dP = -np.power(1+np.exp(2*logP), -2) * np.exp(2*logP)*2

    resultDx[:,0,0] = -theta[1]*expMminusP
    resultDx[:,1,0] = theta[1]*expMminusP
    resultDx[:,2,0] = -theta[0]*np.exp(logH)

    resultDx[:,0,1] = theta[4]*np.multiply(np.exp(-logM),dP)
    resultDx[:,1,1] = -theta[4]*np.exp(-logM)/(1+np.exp(2*logP))

    resultDx[:,0,2] = -theta[0]*np.exp(logP) + theta[5]*np.multiply(np.exp(-logH),dP)
    resultDx[:,2,2] = -theta[5]*np.exp(-logH)/(1+np.exp(2*logP))

    return resultDx


def hes1logmodelDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    logP = (x[:, 0])
    logM = (x[:, 1])
    logH = (x[:, 2])

    resultDtheta[:,0,0] = -np.exp(logH)
    resultDtheta[:,1,0] = np.exp(logM-logP)
    resultDtheta[:,2,0].fill(-1)

    resultDtheta[:,3,1].fill(-1)
    resultDtheta[:,4,1] = np.exp(-logM)/(1+np.exp(2*logP))

    resultDtheta[:,0,2] = -np.exp(logP)
    resultDtheta[:,5,2] = np.exp(-logH)/(1+np.exp(2*logP))
    resultDtheta[:,6,2].fill(-1)

    return resultDtheta


yTest = np.random.rand( np.shape(y)[0], np.shape(y)[1])
thetaTest = np.random.rand( np.shape(true_theta)[0])
testDynamicalModel(hes1logmodelOde, hes1logmodelDx, hes1logmodelDtheta, 
                   "Hes1 log", yTest, thetaTest, tvecObs)

hes1_system = ode_system("Hes1-log-python", hes1logmodelOde, hes1logmodelDx, hes1logmodelDtheta,
                         thetaLowerBound=np.ones(7) * 0, thetaUpperBound=np.ones(7) * np.inf)

control=dict(
    sigma = np.array([0.15, 0.15, np.NaN]),
    useFixedSigma = True
)

yIn = np.append(tvecObs.reshape(-1,1), y, axis = 1)
hes1result = MagiSolver(y=yIn, odeModel=hes1_system, control=control)

# Traceplots of parameters and log-posterior
theta_names = ["a", "b", "c", "d", "e", "f", "g"]
for j in range(hes1result['theta'].shape[0]):
    plt.plot(hes1result['theta'][j,:])
    plt.title(theta_names[j])
    plt.show()
plt.plot(hes1result['lp'])
plt.title('log-post')
plt.show()

# Parameter estimates
theta_est = np.vstack((np.mean(hes1result['theta'], axis=-1),
                   np.quantile(hes1result['theta'], q=[0.025, 0.975], axis=-1)))
pd.DataFrame(theta_est, columns=theta_names, index=['Mean', '2.5%', '97.5%'])


# Inferred trajectories
comp_names = ["P (17 observations)", "M (16 observations)", "H (unobserved)"];
xMean = np.exp(np.mean(hes1result['xsampled'], axis=-1))
xLB = np.exp(np.quantile(hes1result['xsampled'], 0.025, axis=-1))
xUB = np.exp(np.quantile(hes1result['xsampled'], 0.975, axis=-1))
for j in range(xMean.shape[0]):
    plt.plot(tvecObs, xMean[j,:], label='inferred trajectory')
    plt.fill_between(tvecObs, xLB[j,:], xUB[j,:],
                    alpha=0.2, label='95% credible interval')
    plt.scatter(tvecObs, np.exp(y[:,j]), label='noisy observations', c='black')
    plt.plot(tvecOde, sol.y[j, :], label='truth')
    plt.title(comp_names[j])
    plt.xlabel('Time')
    plt.legend()
    plt.show()


## Fitzhugh-Nagumo equations

FNdat = pd.read_csv("../data/FN-sim.csv")

y_I0 = setDiscretization(FNdat.to_numpy(), by = 0.5)

y_I1 = setDiscretization(y_I0, level=1)
y_I2 = setDiscretization(y_I0, level=2)
y_I3 = setDiscretization(y_I0, level=3)


def fnmodelOde(theta, x, tvec):
    V = x[:, 0]
    R = x[:, 1]
    Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
    Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
    result = np.stack([Vdt, Rdt], axis=1)
    return result

def fnmodelDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDx[:,0,0] = theta[2] * (1.0 - np.square(V))
    resultDx[:,1,0] = theta[2]
    resultDx[:,0,1] = -1.0 / theta[2]
    resultDx[:,1,1] = -1.0*theta[1]/theta[2]
    return resultDx

def fnmodelDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
    resultDtheta[:,0,1] = 1.0 / theta[2]
    resultDtheta[:,1,1] = -R / theta[2]
    resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)
    return resultDtheta

fn_system = ode_system("FN-python", fnmodelOde, fnmodelDx, fnmodelDtheta,
                       thetaLowerBound=np.array([0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf]))

FNres0 = MagiSolver(y = y_I0, odeModel=fn_system, control=dict(niterHmc = 10000))
FNres1 = MagiSolver(y = y_I1, odeModel=fn_system, control=dict(niterHmc = 10000))
FNres2 = MagiSolver(y = y_I2, odeModel=fn_system, control=dict(niterHmc = 10000))
FNres3 = MagiSolver(y = y_I3, odeModel=fn_system, control=dict(niterHmc = 10000, nstepsHmc=1000))

def FNpar_table(res):
    allpar = np.vstack((res['theta'], res['sigma']))
    par_est = np.vstack((np.mean(allpar, axis=-1),
                   np.quantile(allpar, q=[0.025, 0.975], axis=-1)))
    return pd.DataFrame(par_est, columns=["a", "b", "c", "sigmaV", "sigmaR"],
                        index=['Mean', '2.5%', '97.5%'])
    
FNpar_table(FNres0)
FNpar_table(FNres1)
FNpar_table(FNres2)
FNpar_table(FNres3)

# Reconstructed trajectories
tvecOde = np.linspace(0,20,2001)
def FNcalcTraj(res):
    x0_est = np.mean(res['xsampled'], axis=-1)[:,0]
    theta_est = np.mean(res['theta'], axis=-1)
    sol = solve_ivp(lambda t, y: fnmodelOde(theta_est, y.transpose(), t).transpose(),
                t_span=[0, 20], y0=x0_est, t_eval=tvecOde, vectorized=True)
    return sol['y']

resList = [FNres0, FNres1, FNres2, FNres3]
FNtr = list(map(FNcalcTraj, resList))
FNtr_labels = ['I0', 'I1', 'I2', 'I3']

plt.scatter(FNdat['time'], FNdat['V'], c='black')
for j in range(len(FNtr)):
    plt.plot(tvecOde, FNtr[j][0,:], label=FNtr_labels[j])
plt.ylabel('V')
plt.xlabel('Time')
plt.legend()
plt.show()

plt.scatter(FNdat['time'], FNdat['R'], c='black')
for j in range(len(FNtr)):
    plt.plot(tvecOde, FNtr[j][1,:], label=FNtr_labels[j])
plt.ylabel('R')
plt.xlabel('Time')
plt.legend()
plt.show()

# RMSDs
FN_rmsd = lambda x: np.sqrt(np.mean((x[:,np.in1d(tvecOde,FNdat['time'])].transpose() - FNdat[['V','R']])**2, axis=0))
pd.DataFrame(list(map(FN_rmsd, FNtr)), index = FNtr_labels).T


## HIV time-dependent example

def hivtdmodelOde(theta, x, tvec):
    TU = x[:,0]
    TI = x[:,1]
    V = x[:,2]

    lambda_val = theta[0]
    rho = theta[1]
    delta = theta[2]
    N = theta[3]
    c = theta[4]

    eta = 9e-5 * (1 - 0.9 * np.cos(np.pi * tvec / 1000))

    result = np.zeros_like(x)
    result[:,0] = lambda_val - rho * TU - eta * TU * V
    result[:,1] = eta * TU * V - delta * TI
    result[:,2] = N * delta * TI - c * V

    return result

def hivtdmodelDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])

    TU = x[:,0]
    TI = x[:,1]
    V = x[:,2]

    lambda_val = theta[0]
    rho = theta[1]
    delta = theta[2]
    N = theta[3]
    c = theta[4]

    eta = 9e-5 * (1 - 0.9 * np.cos(np.pi * tvec / 1000))

    resultDx[:,0,0] = -rho - eta * V
    resultDx[:,1,0] = 0
    resultDx[:,2,0] = -eta * TU

    resultDx[:,0,1] = eta * V
    resultDx[:,1,1] = -delta
    resultDx[:,2,1] = eta * TU

    resultDx[:,0,2] = 0
    resultDx[:,1,2] = N * delta
    resultDx[:,2,2] = -c

    return resultDx

def hivtdmodelDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])

    TU = x[:,0]
    TI = x[:,1]
    V = x[:,2]

    lambda_val = theta[0]
    rho = theta[1]
    delta = theta[2]
    N = theta[3]
    c = theta[4]

    eta = 9e-5 * (1 - 0.9 * np.cos(np.pi * tvec / 1000))

    resultDtheta[:,0,0] = 1
    resultDtheta[:,1,0] = -TU
    resultDtheta[:,2,0] = 0
    resultDtheta[:,3,0] = 0
    resultDtheta[:,4,0] = 0

    resultDtheta[:,0,1] = 0
    resultDtheta[:,1,1] = 0
    resultDtheta[:,2,1] = -TI
    resultDtheta[:,3,1] = 0
    resultDtheta[:,4,1] = 0

    resultDtheta[:,0,2] = 0
    resultDtheta[:,1,2] = 0
    resultDtheta[:,2,2] = N * TI
    resultDtheta[:,3,2] = delta * TI
    resultDtheta[:,4,2] = -V

    return resultDtheta

hiv_time_dependent_system = ode_system("hiv-time-dependent-python", hivtdmodelOde, hivtdmodelDx, hivtdmodelDtheta,
                                       thetaLowerBound=np.array([0,0,0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf, np.inf, np.inf]))

true_theta = [36, 0.108, 0.5, 1000, 3] # lambda, rho, delta, N, c
true_x0 = [600, 30, 1e5] # TU, TI, V initial values
true_sigma = [np.sqrt(10), np.sqrt(10), 10] # noise levels
tvecObs = np.linspace(0, 20, 101) #observation times

sol = solve_ivp(lambda t, y: hivtdmodelOde(true_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecObs[-1]], y0=true_x0, t_eval=tvecObs, vectorized=True)

np.random.seed(12321)

y = sol.y.transpose().copy()

for j in range(np.shape(y)[1]):
    y[:, j] +=  np.random.normal(0, true_sigma[0], np.shape(y)[0])

compnames = ["TU", "TI", "V"]
complabels = ["Concentration", "Concentration", "Load"]

for j in range(np.shape(y)[1]):
    plt.plot(tvecObs, sol.y[j,:])
    plt.scatter(tvecObs, y[:,j], facecolors='none', edgecolors='black')
    plt.title(compnames[j])
    plt.ylabel(complabels[j])
    plt.xlabel('Time')
    plt.show()


# use gpsmoothing to determine phi/sigma
phiEst = np.zeros(shape=[2, np.shape(y)[1]])
sigmaInit = np.zeros(np.shape(y)[1])

for j in range(np.shape(y)[1]):
    hyperparam = gpsmoothing(y[:,j], tvecObs)
    phiEst[:,j] = hyperparam['phi']
    sigmaInit[j] = hyperparam['sigma']

phiEst
sigmaInit

# override phi/sigma for V (3rd) component
phiEst[:,2] = [1e7, 0.5]
sigmaInit[2] = 100.0

yFull = np.append(tvecObs.reshape(-1,1), y, axis = 1)
y_I = setDiscretization(yFull, level=1)

control=dict(
    phi=phiEst,
    sigma=sigmaInit,
)

HIVresult = MagiSolver(y=y_I, odeModel=hiv_time_dependent_system, control=control)

# Parameter estimates
theta_names = ["lambda", "rho", "delta", "N", "c"]
theta_est = np.vstack((np.mean(HIVresult['theta'], axis=-1),
                   np.quantile(HIVresult['theta'], q=[0.025, 0.975], axis=-1)))
pd.DataFrame(theta_est, columns=theta_names, index=['Mean', '2.5%', '97.5%'])

# Inferred trajectories
xMean = np.mean(HIVresult['xsampled'], axis=-1)
xLB = np.quantile(HIVresult['xsampled'], 0.025, axis=-1)
xUB = np.quantile(HIVresult['xsampled'], 0.975, axis=-1)
for j in range(xMean.shape[0]):
    plt.plot(y_I[:,0], xMean[j,:], label='inferred trajectory', linewidth=2)
    plt.fill_between(y_I[:,0], xLB[j,:], xUB[j,:],
                    alpha=0.2, label='95% credible interval')
    plt.plot(tvecObs, sol.y[j, :], label='truth')
    plt.title(compnames[j])
    plt.ylabel(complabels[j])
    plt.xlabel('Time')
    plt.legend()
    plt.show()
