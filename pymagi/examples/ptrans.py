import numpy as np
from arma import ode_system
from magi import MagiSolver


def fOde(theta, x, tvec):
    resultdt = np.zeros(shape=np.shape(x))

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultdt[:, 0] = -theta[0]*S - theta[1] * S * R + theta[2] * RS
    resultdt[:, 1] = theta[0]*S
    resultdt[:, 2] = -theta[1]*S*R + theta[2]*RS + theta[4] * RPP / (theta[5]+RPP)
    resultdt[:, 3] = theta[1]*S*R - theta[2]* RS - theta[3]*RS
    resultdt[:, 4] = theta[3]*RS - theta[4] * RPP / (theta[5]+RPP)

    return resultdt


def fOdeDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultDx[:, 0, 0] = -theta[0] - theta[1] * R
    resultDx[:, 2, 0] = -theta[1] * S
    resultDx[:, 3, 0].fill(theta[2])

    resultDx[:, 0, 1].fill(theta[0])

    resultDx[:, 0, 2] = -theta[1]*R
    resultDx[:, 2, 2] = -theta[1]*S
    resultDx[:, 3, 2].fill(theta[2])
    resultDx[:, 4, 2] =  theta[4] * theta[5] / np.square(theta[5] + RPP)

    resultDx[:, 0, 3] = theta[1]*R
    resultDx[:, 2, 3] = theta[1]*S
    resultDx[:, 3, 3].fill(-theta[2] - theta[3])

    resultDx[:, 3, 4].fill(theta[3])
    resultDx[:, 4, 4] = -theta[4] * theta[5] /  np.square(theta[5] + RPP)

    return resultDx


def fOdeDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])

    S = x[:, 0]
    dS = x[:, 1]
    R = x[:, 2]
    RS = x[:, 3]
    RPP = x[:, 4]

    resultDtheta[:, 0, 0] = -S
    resultDtheta[:, 1, 0] = -S*R
    resultDtheta[:, 2, 0] = RS

    resultDtheta[:, 0, 1] = S

    resultDtheta[:, 1, 2] = -S*R
    resultDtheta[:, 2, 2] = RS
    resultDtheta[:, 4, 2] = RPP / (theta[5]+RPP)
    resultDtheta[:, 5, 2] = -theta[4] * RPP / np.square(theta[5]+RPP)

    resultDtheta[:, 1, 3] = S*R
    resultDtheta[:, 2, 3] = -RS
    resultDtheta[:, 3, 3] = -RS

    resultDtheta[:, 3, 4] = RS
    resultDtheta[:, 4, 4] = - RPP / (theta[5]+RPP)
    resultDtheta[:, 5, 4] = theta[4] * RPP / np.square(theta[5]+RPP)

    return resultDtheta


ptrans_system = ode_system("PTrans-python", fOde, fOdeDx, fOdeDtheta,
                           thetaLowerBound=np.ones(6) * 0, thetaUpperBound=np.ones(6) * 4)


ydataTruth = [[1, 0.588261834720057, 0.405587021811379,
               0.233954596382738, 0.185824926227245, 0.121529475508475, 0.0660579216704765,
               0.0232239721559163, 0.00753621476608807, 0.000635757067732186,
               4.4828522151875e-05, 2.92691291637857e-06, 1.85430809432099e-07,
               7.28853967992039e-10, 2.90513174227738e-12],
              [0, 0.053266895650711,
               0.0873622910225387, 0.130427267370046, 0.145032917209717, 0.166173447332274,
               0.185270502887831, 0.199691529407793, 0.204604196852704, 0.20659618691378,
               0.206753576566759, 0.206764363427542, 0.206765059920321, 0.206765106622966,
               0.206765106806669],
              [1, 0.642586847997489, 0.498289607509476,
               0.384851880112798, 0.360672689559933, 0.337963962897698, 0.334437371299282,
               0.362606647434368, 0.408318304747127, 0.512250740799807, 0.61245271751103,
               0.702776887221291, 0.78106230356887, 0.896447938708228, 0.958939507477765
               ],
              [0, 0.301777886330572, 0.349662193053065, 0.28406917802038,
               0.239159189174826, 0.162847399043611, 0.0890984548705512, 0.0329795416265298,
               0.0122844593001908, 0.00151121723113409, 0.000149977389483994,
               1.26910389636527e-05, 9.71682989611335e-07, 4.82588798220601e-09,
               2.14807760018722e-11],
              [0, 0.0556352656719387, 0.152048199437459,
               0.331078941866822, 0.400168121265241, 0.499188638058692, 0.576464173830167,
               0.604413810939102, 0.579397235952683, 0.48623804196906, 0.387397305099487,
               0.297210421739746, 0.218936724748142, 0.103552056465885, 0.0410604925007539
               ]]
ydataTruth = np.array(ydataTruth).transpose()

ydata = [[1.00003, 0.5887, 0.40509, 0.23316, 0.18567,
          0.12149, 0.06848, 0.024, 0.00674, 0.00109, 0.00044, -0.00099,
          -0.00039, 0.00012, 0.00061],
         [0.00046, 0.05223, 0.08721,
          0.1315, 0.14517, 0.16616, 0.18425, 0.19958, 0.20584, 0.20627,
          0.2076, 0.20519, 0.20674, 0.20635, 0.20808],
         [0.99956,
          0.64249, 0.49785, 0.38488, 0.36059, 0.3376, 0.33559, 0.36219,
          0.40655, 0.51279, 0.61291, 0.70233, 0.77953, 0.89698, 0.95844
          ],
         [0.00244, 0.3011, 0.35024, 0.2845, 0.23995, 0.16463,
          0.08887, 0.03159, 0.01217, -0.00032, -0.00066, -0.00087, 0.00052,
          0.00104, 0.00054],
         [-0.00172, 0.05523, 0.15238, 0.33098,
          0.40037, 0.49987, 0.5774, 0.60517, 0.57964, 0.48722, 0.38806,
          0.29898, 0.21757, 0.1039, 0.04104]]
ydata = np.array(ydata).transpose()
tvecObs = [0, 1, 2, 4, 5, 7, 10, 15, 20, 30, 40, 50, 60, 80, 100]
yFillI0 = np.ndarray([101, 5])
yFillI0.fill(np.nan)
tvecI0 = np.linspace(0, 100, num=101)
for j in range(5):
    yFillI0[:, j] = np.interp(tvecI0, tvecObs, ydata[:, j])


yFull = np.ndarray([201, 5])
yFull.fill(np.nan)
tvecFull = np.linspace(0, 100, num=201)
yFull[[x in tvecObs for x in tvecFull], :] = ydata

xInitExogenous = np.zeros_like(yFull)
for j in range(5):
    xInitExogenous[:, j] = np.interp(tvecFull, tvecObs, ydata[:, j])


### Optimize phi first using equally spaced intervals of 1, i.e., 0,1...,100.

control = dict(
    nstepsHmc = 100,
    niterHmc = 2,
    bandSize = 40,
    burninRatio = 0,
)

hyperInit = MagiSolver(y=yFillI0, odeModel=ptrans_system, tvec=tvecI0, control=control)
sigmaUsed = hyperInit['sigma'][:, 0]

# NEW (Aug 11) ----- plug in sigma estimate and re-estimate phi
control = dict(
    nstepsHmc = 100,
    niterHmc = 2,
    bandSize = 40,
    burninRatio = 0,
    sigma = sigmaUsed,
)
hyperInit = MagiSolver(y=yFillI0, odeModel=ptrans_system, tvec=tvecI0, control=control)
phiUsed = hyperInit['phi']

# sampling
control = dict(
    nstepsHmc = 100,
    niterHmc = 20001,
    bandSize = 40,
    burninRatio = 0.5,
    sigma = sigmaUsed,
    phi = phiUsed,
    xInit = xInitExogenous,
)

result = MagiSolver(y=yFull, odeModel=ptrans_system, tvec=tvecFull, control=control)


inferred_trajectory = np.mean(result['xsampled'], axis=-1)
inferred_theta = np.mean(result['theta'], axis=-1)
