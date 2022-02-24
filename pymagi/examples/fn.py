import numpy as np
from arma import ode_system
from magi import MagiSolver
from scipy.integrate import solve_ivp


def fOde(theta, x, tvec):
    V = x[:, 0]
    R = x[:, 1]
    Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
    Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
    result = np.stack([Vdt, Rdt], axis=1)
    return result

def fOdeDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDx[:,0,0] = theta[2] * (1 - np.square(V))
    resultDx[:,1,0] = theta[2]
    resultDx[:,0,1] = -1.0 / theta[2]
    resultDx[:,1,1] = -1.0*theta[1]/theta[2]
    return resultDx

def fOdeDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
    resultDtheta[:,0,1] = 1.0 / theta[2]
    resultDtheta[:,1,1] = -R / theta[2]
    resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)
    return resultDtheta

fn_system = ode_system("FN-python", fOde, fOdeDx, fOdeDtheta,
                       thetaLowerBound=np.array([0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf]))

true_theta = [0.2,0.2,3]
true_x0 = [-1, 1]
true_sigma = [0.2, 0.2]

tvecObs = np.linspace(0, 20, num=41)
tvecFull = np.linspace(0, 20, num=161)

sol = solve_ivp(lambda t, y: fOde(true_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecFull[-1]], y0=true_x0, t_eval=tvecObs, vectorized=True)

ydataTruth = sol.y

ydataTruth = np.array(ydataTruth).transpose()

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

SEED = np.random.randint(1, 100000)
np.random.seed(SEED)
ydataV = ydataTruth[:, 0] + np.random.normal(0, true_sigma[0], ydataTruth[:, 0].size)
ydataR = ydataTruth[:, 1] + np.random.normal(0, true_sigma[1], ydataTruth[:, 1].size)

ydata = np.stack([np.array(ydataV), np.array(ydataR)], axis=1)
yFull = np.ndarray([161, 2])
yFull.fill(np.nan)
yFull[np.linspace(0, 160, num=41).astype(int), :] = ydata

xInitExogenous = np.zeros_like(yFull)
for j in range(2):
    xInitExogenous[:, j] = np.interp(tvecFull, tvecObs, ydata[:, j])

control=dict(
    nstepsHmc = 100,
    niterHmc = 20001,
    stepSizeFactor = 0.06,
    xInit = xInitExogenous,
)

result = MagiSolver(y=yFull, odeModel=fn_system, tvec=tvecFull, control=control)

inferred_trajectory = np.mean(result['xsampled'], axis=-1)
inferred_theta = np.mean(result['theta'], axis=-1)
np.savetxt("fn_inferred_theta_seed{}.csv".format(SEED), inferred_theta)
np.savetxt("fn_inferred_trajectory_seed{}.csv".format(SEED), inferred_trajectory)
np.savetxt("fn_inferred_sigma_seed{}.csv".format(SEED), np.mean(result['sigma'], axis=-1))

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
