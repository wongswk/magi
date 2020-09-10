addpath('models/');

% Set location of local libcmagi.so before starting MATLAB
% (if not installed system-wide)
% e.g., export LD_LIBRARY_PATH="../cmagi/"

config.nobs = 33;
config.noise = [0.15 0.15 0.1];
config.kernel = "generalMatern";
config.seed = rand(1)*1e7;
config.loglikflag = "withmeanBand";
config.bandsize = 20;
config.hmcSteps = 500;
config.n_iter = 20001;
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

config.priorTemperature = 3;

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

% linear interpolation for xInit
xInit = xsim(:,1)';
xInit = horzcat(xInit',interp1(xsim_obs(1:2:config.nobs,1),xsim_obs(1:2:config.nobs,2),xInit,'linear')',interp1(xsim_obs(2:2:config.nobs,1),xsim_obs(2:2:config.nobs,3),xInit,'linear', 'extrap')',xsim_obs(:,4));

% cpp inference ----------------------------
hes1model.fOde = @hes1logmodelODE;
hes1model.fOdeDx = @hes1logmodelDx;
hes1model.fOdeDtheta= @hes1logmodelDtheta;
hes1model.thetaLowerBound= [0 0 0 0 0 0 0];
hes1model.thetaUpperBound= [Inf Inf Inf Inf Inf Inf Inf];

[samplesCpp, phiUsed] = solveMagi( xsim(:,2:size(xsim,2)), hes1model, xsim(:,1)', pram_true.sigma, [], [], [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

burnin = round(config.n_iter*config.burninRatio);

outCpp = samplesCpp;

gpode.lliklist = outCpp(1,(burnin+1):config.n_iter)';
gpode.xOut = outCpp(2:( 1+(size(xsim,2)-1)*size(xsim,1)),(burnin+1):config.n_iter)';
xLen = (size(xsim,2)-1)*size(xsim,1);
gpode.thOut = outCpp((2+(size(xsim,2)-1)*size(xsim,1)):((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)-1),(burnin+1):config.n_iter)';
gpode.sigmaOut = outCpp(((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)):size(outCpp,1),(burnin+1):config.n_iter)';

% Inferred trajectory and parameter estimates
xEst = reshape(mean(gpode.xOut), [config.nobs 3]);
thetaEst = mean(gpode.thOut);

% Sampled trajectories
for i=1:(size(xsim,2)-1)
    qlim = zeros(2, xLen/(size(xsim,2)-1));
    k=0;
    for j=(1+(i-1)*xLen/(size(xsim,2)-1) ):(xLen/(size(xsim,2)-1)*i)
        k = k+1;
        qlim(:,k) = quantile(gpode.xOut(:,j),[0.025 0.975]);
    end

    subplot(1,size(xsim,2)-1,i)
    fill([xsim(:,1)' fliplr(xsim(:,1)')],[exp(qlim(1,:)) exp(fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(xsim(:,1)', exp(xsim(:,i+1)),'Marker','*')
    hold off;
end

% Histograms of parameters
for i=1:size(gpode.thOut,2)
    subplot(1, size(gpode.thOut,2), i);
    histogram(gpode.thOut(:,i));
    xline(pram_true.theta(i), 'LineWidth', 2, 'Color', 'red');
end
    
% Look at whether these estimates are reasonable for reconstructing
% trajectories using ODE solver
[foo, xRecons]=ode45(@hes1logmodelODEsolve,times,xEst(1,:),[],thetaEst);
plot(times,xRecons);
hold on;
plot(xsim(:,1)', xsim_obs(:,2),'Marker','o','MarkerEdgeColor','blue');
plot(xsim(:,1)', xsim_obs(:,3),'Marker','o','MarkerEdgeColor','red');
hold off;


outFileName =  sprintf('%s-seed%d-noise%.3f', config.modelName, config.seed, config.noise(1));
save( sprintf('%s.mat', outFileName));

