addpath('models/');

% Set location of local libcgpds.so before starting MATLAB 
% (if not installed system-wide)
% e.g., export LD_LIBRARY_PATH="../cmagi/"

config.nobs = 15;
config.noise = 0.01 * ones(1,5); % 0.01 high noise, 0.001 low noise
config.t_end = 100;
config.linfillspace = 0.5;
config.kernel = "generalMatern";
config.seed = rand(1)*1e7;
config.loglikflag = "withmeanBand";
config.bandsize = 40;
config.hmcSteps = 100;
config.n_iter = 20001;
config.burninRatio = 0.50;
config.stepSizeFactor = 0.01;
config.modelName = "PTrans";
config.max_epoch = 1;
config.linearizexInit = true;
config.useMean = true;
config.useBand = true;
config.useFrequencyBasedPrior = true;
config.useScalerSigma = false;
config.useFixedSigma = false;

config.ndis = config.t_end / config.linfillspace + 1;
config.priorTemperature = config.ndis / config.nobs;

config.priorTemperature

pram_true.theta = [0.07, 0.6,0.05,0.3,0.017,0.3];
pram_true.x0 = [1,0,1,0,0];
pram_true.sigma=config.noise;

times = 0:0.1:config.t_end;
[foo, xtrue]=ode45(@ptransmodelODEsolve,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 

xsim = [0,1,2,4,5,7,10,15,20,30,40,50,60,80,100];
xsim = horzcat(xsim',interp1(times,xtrue(:,2),xsim)',interp1(times,xtrue(:,3),xsim)',interp1(times,xtrue(:,4),xsim)',interp1(times,xtrue(:,5),xsim)',interp1(times,xtrue(:,6),xsim)');

rng(config.seed);

for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

xsim_obs = xsim( linspace(1, size(xsim,1), config.nobs),:);

fillC = 0:config.linfillspace:config.t_end;
xsim = zeros(length(fillC),6);
xsim(:,:) = NaN;
xsim(:,1) = fillC';
for i=1:size(xsim,1)
    [isin, loc] = ismember(xsim(i,1), xsim_obs(:,1));
    if isin
        xsim(i,2:6) = xsim_obs(loc,2:6);
    end
end

if config.linearizexInit
    xsimt = xsim(:,1)';
    xsimt = horzcat(xsimt',interp1(xsim_obs(:,1),xsim_obs(:,2),xsimt)',interp1(xsim_obs(:,1),xsim_obs(:,3),xsimt)',interp1(xsim_obs(:,1),xsim_obs(:,4),xsimt)',interp1(xsim_obs(:,1),xsim_obs(:,5),xsimt)',interp1(xsim_obs(:,1),xsim_obs(:,6),xsimt)');
    exoxInit = xsimt(:,2:size(xsim,2));
else
    exoxInit = [];
end

% cpp inference ----------------------------
ptmodel.fOde = @ptransmodelODE;
ptmodel.fOdeDx = @ptransmodelDx;
ptmodel.fOdeDtheta= @ptransmodelDtheta;
ptmodel.thetaLowerBound= [0 0 0 0 0 0];
ptmodel.thetaUpperBound= [4 4 4 4 4 4];

% get initial phi and sigma estimate
[samplesCpp, phiUsed] = solveGpds( exoxInit(ismember(xsim(:,1),0:1:100),:), ptmodel, 0:1:100, [], [], [], [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, 2, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

outCpp = samplesCpp;
outCpp(:,1) = [];
sigmaUsed = outCpp((size(outCpp,1)-size(xsim,2)+2):size(outCpp,1),:)';

% optimize phi
[samplesCpp, phiUsed] = solveGpds( exoxInit(ismember(xsim(:,1),0:1:100),:), ptmodel, 0:1:100, sigmaUsed, [], [], [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, 2, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

% HMC sampling
[samplesCpp, phiUsed] = solveGpds( xsim(:,2:size(xsim,2)), ptmodel, xsim(:,1)', sigmaUsed, phiUsed, exoxInit, [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

burnin = round(config.n_iter*config.burninRatio);

outCpp = samplesCpp;

gpode.lliklist = outCpp(1,(burnin+1):config.n_iter)';
gpode.xOut = outCpp(2:( 1+(size(xsim,2)-1)*size(xsim,1)),(burnin+1):config.n_iter)';
xLen = (size(xsim,2)-1)*size(xsim,1);
gpode.thOut = outCpp((2+(size(xsim,2)-1)*size(xsim,1)):((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)-1),(burnin+1):config.n_iter)';
gpode.sigmaOut = outCpp(((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)):size(outCpp,1),(burnin+1):config.n_iter)';

% Inferred trajectory and parameter estimates
xEst = reshape(mean(gpode.xOut), [config.ndis size(xsim,2)-1]);
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
    fill([xsim(:,1)' fliplr(xsim(:,1)')],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(xsim(:,1)', (xsim(:,i+1)),'Marker','*')
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
[foo, xRecons]=ode45(@ptransmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:(size(xsim,2)-1)
    subplot(1,size(xsim,2)-1,i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end


outFileName =  sprintf('%s-seed%d-noise%.3f', config.modelName, config.seed, config.noise(1));
save( sprintf('%s.mat', outFileName));

