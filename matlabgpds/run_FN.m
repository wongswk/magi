addpath('models/');

config.nobs = 41;
config.noise = [0.2, 0.2];
config.kernel = "generalMatern";
config.seed = 1365546660; %(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
config.loglikflag = "withmeanBand";
config.bandsize = 20;
config.hmcSteps = 500;
config.n_iter = 20001;
config.burninRatio = 0.50;
config.stepSizeFactor = 0.06;
config.filllevel = 2;
config.t_end = 20;
config.modelName = "FN";
config.max_epoch = 1;
config.useMean = true;
config.useBand = true;
config.useFrequencyBasedPrior = true;
config.useScalerSigma = false;
config.useFixedSigma = false;

config.ndis = (config.nobs-1)*2^config.filllevel+1;
config.priorTemperature = config.ndis / config.nobs;

pram_true.theta = [0.2,0.2,3];
pram_true.x0 = [-1, 1];
pram_true.sigma=config.noise;

times = 0:0.1:config.t_end;
[foo, xtrue]=ode45(@fnmodelODEsolve,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 

xsim  = linspace(0,config.t_end,config.nobs);
xsim = horzcat(xsim',interp1(times,xtrue(:,2),xsim)',interp1(times,xtrue(:,3),xsim)');


rng(config.seed);

for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

%xsim_obs = xsim( linspace(1, size(xsim,1), config.nobs),:);
xsim_obs = xsim;

xsim = insertNaN(xsim,config.filllevel);

% linear interpolation for xInit
xInit = xsim(:,1)';
xInit = horzcat(interp1(xsim_obs(:,1),xsim_obs(:,2),xInit,'linear')',interp1(xsim_obs(:,1),xsim_obs(:,3),xInit,'linear')');


% cpp inference ----------------------------
fnmodel.fOde = @fnmodelODE;
fnmodel.fOdeDx = @fnmodelDx;
fnmodel.fOdeDtheta= @fnmodelDtheta;
fnmodel.thetaLowerBound= [0 0 0];
fnmodel.thetaUpperBound= [Inf Inf Inf];

[samplesCpp, phiUsed] = solveGpds( xsim(:,2:size(xsim,2)), fnmodel, xsim(:,1)', [], [], xInit, [], [], [], config.priorTemperature, config.priorTemperature, 1, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
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
[foo, xRecons]=ode45(@fnmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:(size(xsim,2)-1)
    subplot(1,size(xsim,2)-1,i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end


outFileName =  sprintf('%s-seed%d-noise%.3f', config.modelName, config.seed, config.noise(1));
save( sprintf('%s.mat', outFileName));


