addpath('models/');
%addpath('export_fig-master/');
%mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/solveGpds.cpp

config.nobs = 33;
config.noise = [0.15 0.15 0.1];
config.kernel = "generalMatern";
config.seed = 1365546660; %(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
config.loglikflag = "withmeanBand";
config.npostplot = 50;
config.bandsize = 20;
config.hmcSteps = 500;
config.n_iter = 20000;
config.burninRatio = 0.50;
config.stepSizeFactor = 0.01;
config.filllevel = 0;
config.modelName = "Hes1-log";
config.async = true;
config.max_epoch = 1;
config.useMean = true;
config.useBand = true;
config.useFrequencyBasedPrior = true;
config.useScalerSigma = false;
config.useFixedSigma =  true;

config.ndis = (config.nobs-1)*2^config.filllevel+1;
config.priorTemperature = config.ndis / config.nobs;

config.useBand = false; % need to recheck this later

pram_true.theta = [0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3];
pram_true.x0 = [log(1.438575), log(2.037488), log(17.90385)];
pram_true.sigma=config.noise;

times = 0:0.01:240;
[foo, xtrue]=ode45(@hes1logmodelODEsolve,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 

rng(config.seed);
xsim = xtrue;
for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

xsim_obs = xsim( linspace(1, size(xsim,1), config.nobs),:);

xsim_obs(:,4) = NaN;
if config.async
  xsim_obs(2:2:config.nobs,2) = NaN;
  xsim_obs(1:2:config.nobs,3) = NaN;
end

xsim_obs = insertNaN(xsim_obs,config.filllevel);
xsim = xsim_obs;

% cpp inference ----------------------------
hes1model.fOde = @hes1logmodelODE;
hes1model.fOdeDx = @hes1logmodelDx;
hes1model.fOdeDtheta= @hes1logmodelDtheta;
hes1model.thetaLowerBound= [0 0 0 0 0 0 0];
hes1model.thetaUpperBound= [Inf Inf Inf Inf Inf Inf Inf];

[samplesCpp, phiUsed] = solveGpds( xsim(:,2:size(xsim,2)), hes1model, xsim(:,1)', pram_true.sigma, [], [], [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

burnin = round(config.n_iter*config.burninRatio);

outCpp = samplesCpp;

gpode.lliklist = outCpp(1,(burnin+1):config.n_iter)';
gpode.xOut = outCpp(2:( 1+(size(xsim,2)-1)*size(xsim,1)),(burnin+1):config.n_iter)';
xLen = (size(xsim,2)-1)*size(xsim,1);
gpode.thOut = outCpp((2+(size(xsim,2)-1)*size(xsim,1)):((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)-1),(burnin+1):config.n_iter)';
gpode.sigmaOut = outCpp(((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)):size(outCpp,1),(burnin+1):config.n_iter)';

%%% Good up to here.

%%% tempered with warm start
xInit = reshape(mean(gpode.xOut), [config.nobs 3]);
thetaInit = mean(gpode.thOut);
phiNoTemperOptimized = phiUsed;

% Look at whether these estimates are reasonable for reproducing
% trajectories using ODE solver
[foo, xEst]=ode45(@hes1logmodelODEsolve,times,xInit(1,:),[],thetaInit);
plot(times,xEst);
hold on;
plot(xsim(:,1)', xsim_obs(:,2),'Marker','o','MarkerEdgeColor','blue');
plot(xsim(:,1)', xsim_obs(:,3),'Marker','o','MarkerEdgeColor','red');
hold off;

[samplesCpp, phiUsed] = solveGpds( xsim(:,2:size(xsim,2)), hes1model, xsim(:,1)', pram_true.sigma, phiNoTemperOptimized, xInit, thetaInit', [], [], config.priorTemperature, config.priorTemperature, 1/16, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

outCpp = samplesCpp;

gpode.lliklist = outCpp(1,(burnin+1):config.n_iter)';
gpode.xOut = outCpp(2:( 1+(size(xsim,2)-1)*size(xsim,1)),(burnin+1):config.n_iter)';
xLen = (size(xsim,2)-1)*size(xsim,1);
gpode.thOut = outCpp((2+(size(xsim,2)-1)*size(xsim,1)):((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)-1),(burnin+1):config.n_iter)';
gpode.sigmaOut = outCpp(((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)):size(outCpp,1),(burnin+1):config.n_iter)';

% final point estimates for theta
thetaHat = mean(gpode.thOut);


