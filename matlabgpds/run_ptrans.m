addpath('models/');
mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/solveGpds.cpp

config.nobs = 15;
config.noise = 0.001 * ones(1,6); % low noise
config.kernel = "generalMatern";
config.seed = 1365546660; %(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
config.loglikflag = "withmeanBand";
config.bandsize = 20;
config.hmcSteps = 100;
config.n_iter = 1001;
config.burninRatio = 0.50;
config.stepSizeFactor = 0.06;
config.filllevel = 2;
config.modelName = "PTrans";
config.temperPrior = true;
config.useFrequencyBasedPrior = true;
config.useScalerSigma = false;
config.useFixedSigma =  false;
config.max_epoch = 10;

config.ndis = (config.nobs-1)*2^config.filllevel+1;
if (config.temperPrior)
  config.priorTemperature = config.ndis / config.nobs;
else
  config.priorTemperature = 1;
end

if config.loglikflag == "withmeanBand"
  config.useMean = true;
  config.useBand = true;
elseif config.loglikflag == "band"
  config.useMean = false;
  config.useBand = true;
elseif config.loglikflag == "withmean"
  config.useMean = true;
  config.useBand = false;
elseif config.loglikflag == "usual"
  config.useMean = false;
  config.useBand = false;
end

pram_true.theta = [0.07, 0.6,0.05,0.3,0.017,0.3];
pram_true.x0 = [1,0,1,0,0];
pram_true.sigma=config.noise;

times = 0:0.1:100;
[foo, xtrue]=ode45(@ptransmodelODEsolve,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 
plot(xtrue(:,1),xtrue(:,2:6));

xsim = [0,1,2,4,5,7,10,15,20,30,40,50,60,80,100];
%xsim  = linspace(0,100,config.nobs);
xsim = horzcat(xsim',interp1(times,xtrue(:,2),xsim)',interp1(times,xtrue(:,3),xsim)',interp1(times,xtrue(:,4),xsim)',interp1(times,xtrue(:,5),xsim)',interp1(times,xtrue(:,6),xsim)');

rng(config.seed);

for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

xsim_obs = xsim( linspace(1, size(xsim,1), config.nobs),:);

hold on;
plot(xsim_obs(:,1), xsim_obs(:, 2:size(xsim,2)), 'o');
hold off;

plot(xsim_obs(:,1), xsim_obs(:, 2:size(xsim,2)), 'o');

xsim = insertNaN(xsim_obs,config.filllevel);

% cpp inference ----------------------------
ptmodel.fOde = @ptransmodelODE;
ptmodel.fOdeDx = @ptransmodelDx;
ptmodel.fOdeDtheta= @ptransmodelDtheta;
ptmodel.thetaLowerBound= [0 0 0 0 0 0];
ptmodel.thetaUpperBound= [Inf Inf Inf Inf Inf Inf];

samplesCpp = solveGpds( xsim(:,2:size(xsim,2)), ptmodel, xsim(:,1)', [], config.priorTemperature, config.priorTemperature, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

