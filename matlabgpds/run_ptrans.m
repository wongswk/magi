addpath('models/');
addpath('export_fig-master/');
%mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/solveGpds.cpp

config.nobs = 15;
config.noise = 0.01 * ones(1,5); % 0.01 high noise, 0.001 low noise
config.t_end = 100;
config.kernel = "generalMatern";
config.seed = 1365546660; %(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
config.loglikflag = "withmeanBand";
config.bandsize = 20;
config.hmcSteps = 100;
config.n_iter = 10001;
config.burninRatio = 0.50;
config.stepSizeFactor = 0.06;
config.filllevel = 3;
config.modelName = "PTrans";
config.temperPrior = true;
config.useFrequencyBasedPrior = true;
config.useScalerSigma = false;
config.useFixedSigma =  false;
config.max_epoch = 10;
config.useExoSigma = true;
config.linearizexInit = true;

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

if config.useExoSigma
    exoSigma = config.noise;
else
    exoSigma = [];
end


xsim = insertNaN(xsim_obs,config.filllevel);

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
ptmodel.thetaUpperBound= [Inf Inf Inf Inf Inf Inf];

samplesCpp = solveGpds( xsim(:,2:size(xsim,2)), ptmodel, xsim(:,1)', exoSigma, [], exoxInit, [], [], config.priorTemperature, config.priorTemperature, char(config.kernel), config.hmcSteps, config.burninRatio, config.n_iter, ...
    config.stepSizeFactor, config.max_epoch, config.bandsize, config.useFrequencyBasedPrior, config.useBand, config.useMean, config.useScalerSigma, config.useFixedSigma, true);

outCpp = samplesCpp;
outCpp(:,1) = [];
lliklist = outCpp(1,:);
xOut = outCpp(2:( 1+(size(xsim,2)-1)*size(xsim,1)),:)';
xLen = (size(xsim,2)-1)*size(xsim,1);
thOut = outCpp((2+(size(xsim,2)-1)*size(xsim,1)):((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)-1),:)';
sigmaOut = outCpp(((2+(size(xsim,2)-1)*size(xsim,1))+length(pram_true.theta)):size(outCpp,1),:)';

outFileName =  sprintf('result-%s-noise%.3f-fill%d-useExoSigma%d-useFixedSigma%d-linearizeInit%d', config.modelName, config.noise(1), config.filllevel, config.useExoSigma, config.useFixedSigma, config.linearizexInit);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 1 1])

fig = figure;
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 11 8.5])

for i=1:(size(xsim,2)-1)
    qlim = zeros(2, xLen/(size(xsim,2)-1));
    k=0;
    for j=(1+(i-1)*xLen/(size(xsim,2)-1) ):(xLen/(size(xsim,2)-1)*i)
        k = k+1;
        qlim(:,k) = quantile(xOut( ((config.n_iter-1) * config.burninRatio):(config.n_iter-1),j),[0.025 0.975]);
    end

    subplot(2,3,i)
    fill([xsim(:,1)' fliplr(xsim(:,1)')],[qlim(1,:) fliplr(qlim(2,:))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(xsim(:,1)',xsim(:,i+1),'Marker','*')
    hold off;
end

subplot(2,3,6)
plot(lliklist)
export_fig(outFileName, '-pdf');

figure;
xsim_len = length(xsim);
xInitr = reshape(xOut(1,:),[xsim_len 5])';
for i=1:5
    subplot(2,3,i)
    plot(xsim_obs(:,1), xsim_obs(:, i+1), 'o');
    hold on;
    plot(xsim(:,1), xInitr(i,:));
    hold off;
end
export_fig(outFileName, '-pdf', '-append');

figure;
for i=1:6
    subplot(2,3,i);
    plot(thOut(:,i));
    yline(pram_true.theta(i), 'LineWidth', 2, 'Color', 'red');
end
export_fig(outFileName, '-pdf', '-append');


figure;
for i=1:5
    subplot(2,3,i);
    plot(sigmaOut(:,i));
end
export_fig(outFileName, '-pdf', '-append');

xtest = zeros(100, length(times), size(xsim,2)-1);
repIter = round(linspace( (config.n_iter-1) * config.burninRatio,config.n_iter-1,100),0);
for i=1:100
    [foo, xtest(i,:,:)]=ode45(@ptransmodelODEsolve,times,xOut(repIter(i),1:size(xsim,1):xLen),[],thOut(repIter(i),:));
end

[maxllik, MAPIndex] = max(lliklist);
[foo, x_MAP]=ode45(@ptransmodelODEsolve,times,xOut(MAPIndex,1:size(xsim,1):xLen),[],thOut(MAPIndex,:));

ode025 = zeros(length(times), size(xsim,2)-1);
ode975 = zeros(length(times), size(xsim,2)-1);
for i=1:length(times)
    for j = 1:(size(xsim,2)-1)
        foo = quantile(xtest(:,i,j),[0.025 0.975]);
        ode025(i,j) = foo(1);
        ode975(i,j) = foo(2);
    end
end

figure;
for i=1:(size(xsim,2)-1)
    subplot(2,3,i)
    fill([times fliplr(times)],[ode025(:,i)' fliplr(ode975(:,i)')],[.9 .9 .9],'LineStyle','none')
    hold on;
    plot(times, ode975(:,i), 'Color', 'blue');
    plot(times, ode025(:,i), 'Color', 'blue');
    plot(times, x_MAP(:,i), 'Color', 'black', 'LineWidth', 2);
    plot(xsim_obs(:,1)',xsim_obs(:,i+1),'Marker','*', 'Color', 'red');
    hold off;
end
export_fig(outFileName, '-pdf', '-append');

save( sprintf('%s.mat', outFileName))

