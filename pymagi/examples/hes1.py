import numpy as np
from arma import ode_system
from magi import MagiSolver
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp


def fOde(theta, x, tvec):
    P = np.exp(x[:, 0])
    M = np.exp(x[:, 1])
    H = np.exp(x[:, 2])
    PMHdt = np.zeros(shape=np.shape(x))
    PMHdt[:,0] = -theta[0]*H + theta[1]*M/P - theta[2]
    PMHdt[:,1] = -theta[3] + theta[4]/(1+np.square(P))/M
    PMHdt[:,2] = -theta[0]*P + theta[5]/(1+np.square(P))/H - theta[6]
    return PMHdt


def fOdeDx(theta, x, tvec):
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

    return resultDx


def fOdeDtheta(theta, x, tvec):
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

    return resultDtheta


hes1_system = ode_system("Hes1-log-python", fOde, fOdeDx, fOdeDtheta,
                         thetaLowerBound=np.ones(7) * 0, thetaUpperBound=np.ones(7) * np.inf)

true_theta = [0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3]
true_x0 = [np.log(1.438575), np.log(2.037488), np.log(17.90385)]
true_sigma=[0.15, 0.15, np.NaN]

tvecFull = np.linspace(0, 240, num=33)

sol = solve_ivp(lambda t, y: fOde(true_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecFull[-1]], y0=true_x0, t_eval=tvecFull, vectorized=True)

ydataTruth = sol.y
ydataTruth = np.array(ydataTruth).transpose()


ydataP = [0.74, 0.78, 1.86, 1.86, 2.2, 1.93, 1.47, 1.03, 0.36,
          0.88, 1.68, 1.97, 2.15, 1.85, 1.8, 1.47, 0.71]
ydataM = [0.91, 0.82, 0.71, -0.11, 0.08, -0.45, -0.05, 0.2,
          0.88, 1.09, 0.3, 0.35, 0.25, -0.23, -0.51, -0.09]

SEED = np.random.randint(1, 100000)
np.random.seed(SEED)
ydataP = ydataTruth[np.linspace(0, 32, num=17).astype(int), 0] + np.random.normal(0, 0.15, 17)
ydataM = ydataTruth[np.linspace(1, 31, num=16).astype(int), 1] + np.random.normal(0, 0.15, 16)

yFull = np.ndarray([33, 3])
yFull.fill(np.nan)
yFull[np.linspace(0, 32, num=17).astype(int), 0] = ydataP
yFull[np.linspace(1, 31, num=16).astype(int), 1] = ydataM

control=dict(
    nstepsHmc = 500,
    niterHmc = 20001,
    sigma = np.array([0.15, 0.15, np.NaN]),
    useFixedSigma = True,
)

result = MagiSolver(y=yFull, odeModel=hes1_system, tvec=tvecFull, control=control)

inferred_trajectory = np.mean(result['xsampled'], axis=-1)
inferred_theta = np.mean(result['theta'], axis=-1)

np.savetxt("hes1log_inferred_theta_seed{}.csv".format(SEED), inferred_theta)
np.savetxt("hes1log_inferred_trajectory_seed{}.csv".format(SEED), inferred_trajectory)
np.savetxt("hes1log_inferred_sigma_seed{}.csv".format(SEED), np.mean(result['sigma'], axis=-1))

# Inferred trajectories visualization
for j in range(inferred_trajectory.shape[0]):
    plt.plot(tvecFull, inferred_trajectory[j,:])
plt.show()

for j in range(inferred_trajectory.shape[0]):
    plt.plot(tvecFull, np.quantile(result['xsampled'], 0.025, axis=-1)[j,:])
    plt.plot(tvecFull, np.quantile(result['xsampled'], 0.975, axis=-1)[j,:])
plt.show()

# Histogram of parameters
for j in range(result['theta'].shape[0]):
    plt.hist(result['theta'][j,:])
    plt.show()

# Look at whether these estimates are reasonable for reconstructing trajectories using ODE solver
sol = solve_ivp(lambda t, y: fOde(inferred_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecFull[-1]], y0=inferred_trajectory[:, 0], t_eval=tvecFull, vectorized=True)

for j in range(inferred_trajectory.shape[0]):
    plt.plot(sol.t, sol.y[j,:])
    plt.scatter(tvecFull, yFull[:, j])
plt.show()
