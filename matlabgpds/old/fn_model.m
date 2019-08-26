%  Self-contained FN model demonstration using C++ defined OdeSystem
%  (withBand currently needs to be disabled, maybe issue with LAPACK linked?)

mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/phisigllikC.cpp
mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/xthetallikM.cpp
mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/xthetasigmaSample.cpp

config.nobs = 41;
config.noise = [0.15, 0.07] * 2;
config.kernel = "generalMatern";
config.seed = 1365546660; %(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
config.npostplot = 50;
config.loglikflag = "withmeanBand";
config.bandsize = 20;
config.hmcSteps = 500;
config.n_iter = 1e4;
config.burninRatio = 0.50;
config.stepSizeFactor = 1;
config.filllevel = 2;
config.modelName = "FN";
config.startXAtTruth = false;
config.startThetaAtTruth = false;
config.startSigmaAtTruth = false;
config.useGPmean = true;
config.forseTrueMean = false;
config.phase2 = false;
config.temperPrior = true;
config.phase3 = false;
config.max_epoch = 10;not working
config.epoch_method = ["mean", "median", "deSolve", "f_x_bar"];
config.epoch_method = config.epoch_method(1);


config.ndis = (config.nobs-1)*2^config.filllevel+1;
if (config.temperPrior)
  config.priorTemperature = config.ndis / config.nobs;
  config.priorString = 'temperPrior';
else
  config.priorTemperature = 1;
  config.priorString = 'unitHeatPrior';
end

% initialize global parameters, true x, simulated x ----------------------------
% baseDir = '/home/s246wong/dynsys-Matlab/';
% outDir =  baseDir + config.modelName + '-' + config.loglikflag +  '-' + config.kernel ...
%                               + '-nobs' + config.nobs +  '-noise' + num2str(config.noise, '%#5.3f_')...
%                               + '-ndis' + config.ndis + '-' + config.priorString + '/';
% 
% system("mkdir -p " + outDir);

pram_true.theta = [0.2,0.2,3];
pram_true.x0 = [-1, 1];
pram_true.phi= [0.9486433, 3.2682434;
                1.9840824, 1.1185157];
pram_true.sigma=config.noise;

times = 0:(20/240):20;

[foo, xtrue]=ode45(@fnmodelODE,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 
plot(xtrue(:,1),xtrue(:,2:3));

xsim  = linspace(0,20,config.nobs);
xsim = horzcat(xsim',interp1(times,xtrue(:,2),xsim)',interp1(times,xtrue(:,3),xsim)');

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

tvec.full = xsim(:,1);
tvec.nobs = xsim_obs(:,1);

foo1 = @(x,y)x-y;
foo = bsxfun(foo1, tvec.full, tvec.full');
r = abs(foo);
r2 = r.^2;
signr = -sign(foo);

foo = bsxfun(foo1, tvec.nobs, tvec.nobs');
r_nobs = abs(foo);
r2_nobs = r_nobs.^2;
signr_nobs = -sign(foo);

% GPsmoothing: marllik+fftNormalprior for phi-sigma ----------------------------
cursigma = NaN(1,size(xsim,2)-1);
curphi = NaN(2, size(xsim,2)-1);

for j=1:(size(xsim,2)-1)
    [priormeanFactor, priorsdFactor] = getFrequencyBasedPrior(xsim_obs(:,j+1));

    desiredMode = priormeanFactor;
    betaRate = fzero(@(betaRate) gamcdf(1, 1 + desiredMode*betaRate, 1/betaRate)-0.95, [1e-3 1e3]);
    alphaRate = 1 + desiredMode*betaRate;

    mlp = @(pars)marlikprior(pars, xsim_obs, r_nobs, config, j, priormeanFactor, priorsdFactor);
    optim_options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true);
    ostart = [std(xsim_obs(:,1+j))/2 max(xsim_obs(:,1))/2 std(xsim_obs(:,1+j))/2]; 
    [xout,fval,exitflag,output,grad] = fminunc(mlp,ostart,optim_options);
    
    cursigma(j) = xout(3);
    curphi(j,:) = xout(1:2);
end


for j=1:(size(xsim,2)-1)
   curCov(j) =  calCov(curphi(:,j), r, signr, config.bandsize);
end

if(config.useGPmean)
  for j=1:(size(xsim,2)-1)
    ydy = getMeanCurve(xsim_obs(:,1), xsim_obs(:,j+1), xsim(:,1), ...
                        curphi(:,j)', cursigma(j));
                        %kerneltype=config$kernel, deriv = TRUE);
    curCov(j).mu = ydy.y;
    curCov(j).dotmu = ydy.dy;
  end
end

% MCMC starting value ----------------------------------------------------------
yobs = xsim(:,2:size(xsim,2));

xsimInit = xsim;
for j=1:(size(xsim,2)-1)
  xsimInit(:,j+1) = interp1( xsim_obs(:,1), xsim_obs(:,j+1), xsim(:,1), 'pchip');
end
xInit = reshape(xsimInit(:,2:size(xsim,2)),[1 size(xsim,1)*(size(xsim,2)-1)]);

thetaInit = ones(1, length(pram_true.theta));


llp = @(theta) thetaoptim(yobs, curCov, cursigma, [xInit theta], char(config.modelName), false );

    [xout,fval,exitflag,output,grad] = fminunc(llp,thetaInit,optim_options);
    thetaInit = xout;
    sigmaInit = cursigma;
xthetasigmaInit = [ xInit, thetaInit, sigmaInit ] ;
stepLowInit = ones(1, length(xthetasigmaInit)) * 0.000035;
stepLowInit = stepLowInit*config.stepSizeFactor;

xId = 1:length(xInit);
thetaId = (max(xId)+1):(max(xId)+length(thetaInit));
sigmaId = (max(thetaId)+1):(max(thetaId)+length(sigmaInit));

xthetasigamSingleSampler = @(xthetasigma, stepSize)  xthetasigmaSample(yobs, curCov, xthetasigma(sigmaId), xthetasigma([xId, thetaId]), ...
                    stepSize, config.hmcSteps, false, char(config.loglikflag), ...
                    config.priorTemperature, char(config.modelName));

chainSamplesOut = chainSampler(config, xthetasigmaInit, xthetasigamSingleSampler, stepLowInit, true);

%%%%%%%%%%%%%% Generate summary figure.
fig = figure;
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 11 8.5])

j=0;
for i=(max(xId)+1):(max(xId)+length(thetaInit))
   j = j+1;
   subplot(2,3,j)
   histogram(chainSamplesOut.xth( (config.n_iter * config.burninRatio):config.n_iter,i)); 
end

for i=1:(size(xsim,2)-1)
    qlim = zeros(2, length(xInit)/(size(xsim,2)-1));
    k=0;
    for j=(1+(i-1)*length(xInit)/(size(xsim,2)-1) ):(length(xInit)/(size(xsim,2)-1)*i)
        k = k+1;
        qlim(:,k) = quantile(chainSamplesOut.xth( (config.n_iter * config.burninRatio):config.n_iter,j),[0.025 0.975]);
    end

    subplot(2,3,3+i)
    fill([tvec.full' fliplr(tvec.full')],[qlim(1,:) fliplr(qlim(2,:))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(tvec.full',yobs(:,i),'Marker','*')
    hold off;
end

print(fig,'results','-dpdf')



