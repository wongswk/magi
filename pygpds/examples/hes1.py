import numpy as np
from arma import ode_system, solve_gpds


def fOde(theta, x):
    P = np.exp(x[:, 0])
    M = np.exp(x[:, 1])
    H = np.exp(x[:, 2])
    PMHdt = np.zeros(shape=np.shape(x))
    PMHdt[:,0] = -theta[0]*H + theta[1]*M/P - theta[2]
    PMHdt[:,1] = -theta[3] + theta[4]/(1+np.square(P))/M
    PMHdt[:,2] = -theta[0]*P + theta[5]/(1+np.square(P))/H - theta[6]
    return PMHdt


def fOdeDx(theta, x):
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


def fOdeDtheta(theta, x):
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


ydataTruth = [[0.363653040293724, 0.577225047385934, 1.02000277711257,
               1.48433447382444, 1.83746499681482, 2.0497916129288, 2.15351882026633,
               2.18288565908666, 2.15985983799388, 2.09691863056832, 2.00035349182604,
               1.87183797482343, 1.70847770782502, 1.50127538160754, 1.23178627383479,
               0.873645719233323, 0.475982900087656, 0.392228800672109, 0.713495479940593,
               1.18288142322575, 1.61980272155199, 1.92385302748331, 2.09504781703248,
               2.17042341778082, 2.18018516555888, 2.14263260748799, 2.06787050409465,
               1.96053452063022, 1.82082102975053, 1.64403607292539, 1.41843207518635,
               1.12172319590024, 0.733526912720838],
              [0.711717676544985,
               0.97115737201471, 1.00131842943306, 0.901328674077051, 0.747756838900062,
               0.575876929124666, 0.399293976052713, 0.224425679429039, 0.0560591262946247,
               -0.100450143557316, -0.237771149639473, -0.345127880481337, -0.406433204640876,
               -0.397823209415185, -0.284374716892247, -0.0187262283156898,
               0.415947665206159, 0.829087233244039, 1.00236396719822, 0.97729572698794,
               0.85296426777186, 0.690778707082951, 0.516330146219065, 0.339755129468047,
               0.166552527408967, 0.00152352981681982, -0.149459154161395, -0.278079713273549,
               -0.372001594392621, -0.412821659998205, -0.373346799828477, -0.214155602559634,
               0.112445464261916],
              [2.88501577351839, 2.80695422810485,
               2.20313246828105, 1.29208378373985, 0.393259918372852, -0.227647398445347,
               -0.552463328411332, -0.669419725875601, -0.646388892760384, -0.522454269563155,
               -0.318936235134039, -0.0456218413158569, 0.297100070888378, 0.718918428439833,
               1.24183779140676, 1.89194481225458, 2.59263722342394, 2.92891381761992,
               2.65212608225424, 1.91585548235356, 0.969317632090725, 0.146611170217348,
               -0.366323730542059, -0.611131674472768, -0.674872789806283, -0.614440501124868,
               -0.461879564657812, -0.234159461824482, 0.062191976517597, 0.42999485046593,
               0.882873661322661, 1.4467562116298, 2.13585225078904]]
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

tvecFull = np.linspace(0, 240, num=33)


result = solve_gpds(
    yFull,
    hes1_system,
    tvecFull,
    sigmaExogenous = np.array([0.15, 0.15, 0.15]),
    phiExogenous = np.array([[]]),
    xInitExogenous = np.array([[]]),
    thetaInitExogenous = np.array([]),
    muExogenous = np.array([[]]),
    dotmuExogenous = np.array([[]]),
    priorTemperatureLevel = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureDeriv = yFull.size/np.sum(np.isfinite(yFull)),
    priorTemperatureObs = 1.0,
    kernel = "generalMatern",
    nstepsHmc = 500,
    burninRatioHmc = 0.5,
    niterHmc = 20001,
    stepSizeFactorHmc = 0.01,
    nEpoch = 1,
    bandSize = 20,
    useFrequencyBasedPrior = True,
    useBand = True,
    useMean = True,
    useScalerSigma = False,
    useFixedSigma = True,
    verbose = True)

print(result['phiUsed'])
print(result['samplesCpp'])

# verify trajectory RMSE
samplesCpp = result['samplesCpp']
llikId = 0
xId = range(np.max(llikId)+1, np.max(llikId)+yFull.size+1)
thetaId = range(np.max(xId)+1, np.max(xId)+7+1)
sigmaId = range(np.max(thetaId)+1, np.max(thetaId)+yFull.shape[1]+1)

burnin = int(20001*0.5)
xsampled = samplesCpp[xId, (burnin+1):]
xsampled = xsampled.reshape([yFull.shape[1], yFull.shape[0], -1])
thetaSampled = samplesCpp[thetaId, (burnin+1):]
sigmaSampled = samplesCpp[sigmaId, (burnin+1):]

from matplotlib import pyplot as plt
for j in range(yFull.shape[1]):
    plt.plot(tvecFull, xsampled[j, :, -1])  # one single sample plot

inferred_trajectory = np.mean(xsampled, axis=-1)
for j in range(yFull.shape[1]):
    plt.plot(tvecFull, inferred_trajectory[j, :])  # inferred trajectory plot


thetaSampled = samplesCpp[thetaId, (burnin+1):]
inferred_theta = np.mean(thetaSampled, axis=-1)
np.savetxt("hes1log_inferred_theta_seed{}.csv".format(SEED), inferred_theta)
np.savetxt("hes1log_inferred_trajectory_seed{}.csv".format(SEED), inferred_trajectory)
np.savetxt("hes1log_inferred_sigma_seed{}.csv".format(SEED), np.mean(samplesCpp[sigmaId, (burnin+1):], axis=-1))
