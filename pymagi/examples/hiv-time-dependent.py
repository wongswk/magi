import numpy as np
from arma import ode_system
from magi import MagiSolver
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt


def fOde(theta, x, tvec):
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

def fOdeDx(theta, x, tvec):
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

def fOdeDtheta(theta, x, tvec):
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

hiv_time_dependent_system = ode_system("hiv-time-dependent-python", fOde, fOdeDx, fOdeDtheta,
                                       thetaLowerBound=np.array([0,0,0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf, np.inf, np.inf]))

true_theta = [36, 0.108, 0.5, 1000, 3] # lambda, rho, delta, N, c
true_x0 = [600, 30, 1e5] # TU, TI, V
true_sigma = [np.sqrt(10), np.sqrt(10), 10]

tvecObs = np.linspace(0, 20, num=101)
tvecFull = np.linspace(0, 20, num=201)

sol = solve_ivp(lambda t, y: fOde(true_theta, y.transpose(), t).transpose(),
                t_span=[0, tvecFull[-1]], y0=true_x0, t_eval=tvecObs, vectorized=True)

ydataTruth = sol.y
ydataTruth = np.array(ydataTruth).transpose()

ydata = [[605.482, 516.696, 466.816, 437.609, 410.107,
          386.236, 363.562, 351.9, 327.346, 320.566, 301.164, 290.443,
          268.17, 257.382, 247.484, 237.934, 229.424, 221.438, 217.612,
          206.664, 204.704, 196.382, 190.519, 193.683, 185.257, 186.968,
          181.519, 176.241, 170.791, 165.105, 166.867, 168.602, 169.359,
          169.918, 170.889, 166.512, 168.321, 171.22, 167.41, 171.782,
          169.916, 170.47, 170.062, 164.391, 165.992, 167.702, 171.332,
          181.286, 180.093, 179.446, 169.309, 181.446, 178.755, 185.912,
          183.024, 181.192, 186.478, 187.177, 186.718, 186.615, 194.321,
          192.726, 199.17, 200.025, 192.407, 203.954, 201.518, 201.224,
          201.932, 211.742, 202.055, 203.69, 216.964, 211.345, 217.715,
          210.954, 220.29, 221.065, 221.351, 227.62, 223.622, 219.565,
          221.877, 227.99, 227.979, 234.202, 227.661, 231.168, 234.044,
          233.532, 233.119, 241.112, 238.959, 243.17, 241.438, 246.203,
          248.597, 252.271, 246.703, 248.263, 256.479],
         [28.828,
          101.493, 138.154, 153.915, 157.061, 165.231, 165.133, 167.071,
          165.987, 165.605, 164.141, 162.585, 163.712, 161.347, 156.114,
          150.853, 145.825, 145.315, 136.092, 132.37, 133.614, 125.023,
          128.136, 118.319, 112.14, 105.574, 102.913, 101.924, 95.475,
          88.811, 87.937, 84.995, 81.062, 75.273, 73.739, 71.83, 65.756,
          64.204, 58.246, 54.373, 53.604, 50.421, 50.975, 51.848, 51.729,
          43.989, 35.588, 39.843, 39.285, 41.489, 36.038, 35.834, 32.324,
          30.295, 28.331, 30.309, 26.301, 25.13, 27.754, 25.368, 30.553,
          23.74, 22.294, 26.699, 18.23, 20.889, 16.237, 11.566, 22.126,
          12.092, 18.38, 19.159, 17.567, 10.487, 18.129, 14.306, 12.576,
          17.264, 7.828, 13.507, 14.951, 9.082, 11.897, 12.76, 9.653, 9.192,
          9.983, 6.67, 8.767, 7.018, 9.48, 7.75, 2.39, 8.006, 9.067, 12.725,
          15.606, 5.506, 3.398, 5.953, 8.474],
         [99987.475, 60396.44,
          42237.232, 33982.597, 30327.029, 28790.103, 28226.011, 28069.724,
          28052.865, 28043.399, 27964.969, 27789.221, 27536.936, 27184.736,
          26770.611, 26271.002, 25728.799, 25114.137, 24436.717, 23748.905,
          23021.303, 22257.192, 21518.631, 20757.727, 19986.504, 19225.924,
          18460.674, 17723.894, 16996.794, 16296.273, 15608.727, 14934.377,
          14291.511, 13672.176, 13083.674, 12498.776, 11958.922, 11438.808,
          10918.491, 10423.493, 9976.698, 9531.598, 9113.649, 8730.256,
          8342.503, 7968.776, 7644.89, 7286.271, 6984.949, 6690.946, 6434.56,
          6153.899, 5913.206, 5637.866, 5413.823, 5200.372, 4994.445, 4808.642,
          4609.432, 4430.109, 4279.908, 4104.128, 3944.593, 3799.133, 3659.732,
          3528.605, 3412.804, 3286.724, 3162.941, 3086.033, 2979.277, 2888.409,
          2782.166, 2680.629, 2592.124, 2527.072, 2442.519, 2366.37, 2301.341,
          2212.542, 2157.973, 2103.708, 2052.807, 1979.577, 1939.42, 1878.018,
          1821.99, 1788.558, 1748.996, 1693.653, 1637.759, 1592.744, 1582.665,
          1522.955, 1491.6, 1459.352, 1427.036, 1417.456, 1351.105, 1370.643,
          1305.584]]
ydata = np.array(ydata).transpose()

SEED = np.random.randint(1, 100000)
np.random.seed(SEED)

ydata[:, 0] = ydataTruth[:, 0] + np.random.normal(0, true_sigma[0], ydataTruth[:, 0].size)
ydata[:, 1] = ydataTruth[:, 1] + np.random.normal(0, true_sigma[1], ydataTruth[:, 1].size)
ydata[:, 2] = ydataTruth[:, 2] + np.random.normal(0, true_sigma[2], ydataTruth[:, 2].size)

yFull = np.ndarray([201, 3])
yFull.fill(np.nan)
yFull[np.linspace(0, 20, num=101).astype(int), :] = ydata

xInitExogenous = np.zeros_like(yFull)
for j in range(3):
    xInitExogenous[:, j] = np.interp(tvecFull, tvecObs, ydata[:, j])

#' manually override estimated hyper-parameters for component 3 because
#' GP smoothing gives bad result for rapidly decreasing curve
phiExogenous = np.array([37322.66, 11271.63, 14299.71, 4.16, 3.02, 1]).reshape([2, ydata.shape[1]])
sigmaInit = np.array([4, 3, 10])


control=dict(
    niterHmc = 20000,
    stepSizeFactor = 0.06,
    xInit = xInitExogenous,
    phi=phiExogenous,
    sigma=sigmaInit,
)

result = MagiSolver(y=yFull, odeModel=hiv_time_dependent_system, tvec=tvecFull, control=control)

inferred_trajectory = np.mean(result['xsampled'], axis=-1)
inferred_theta = np.mean(result['theta'], axis=-1)
np.savetxt("hiv_time_dependent_inferred_theta_seed{}.csv".format(SEED), inferred_theta)
np.savetxt("hiv_time_dependent_inferred_trajectory_seed{}.csv".format(SEED), inferred_trajectory)
np.savetxt("hiv_time_dependent_inferred_sigma_seed{}.csv".format(SEED), np.mean(result['sigma'], axis=-1))

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
