% MATLAB replication script

% To run this script, ensure that path contains 
% - libcmagi.so (or libcmagi.dll on Windows) 
% - 'mex' compiled .cpp files from src/
% addpath('...'); % add path here if necessary

addpath('models/');

%%% Hes1 example

pram_true.theta = [0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3];
pram_true.x0 = [1.439, 2.037, 17.904];
pram_true.sigma = [0.15 0.15 NaN];
nobs = 33;

times = 0:0.01:240;
[foo, xtrue]=ode45(@hes1modelODEsolve,times,pram_true.x0,[],pram_true.theta);

rng(12321);

y = horzcat(foo,xtrue);
y = y( linspace(1, size(y,1), nobs),:);
for j=1:(size(y,2)-1)
  y(:,1+j) = y(:, 1+j) .* exp(normrnd(zeros(size(y,1),1), pram_true.sigma(j)));
end

% Set asynchronous observation schedule
y(:,4) = NaN;
y(2:2:nobs,2) = NaN;
y(1:2:nobs,3) = NaN;

plot(times, xtrue);
legend('P', 'M', 'H', 'AutoUpdate', 'off');
hold on; plot(y(:,1), y(:,2:4), 'Marker','.', 'Color', 'black');
xlabel('Time (min)');
ylabel('Level');
hold off;

% Use log transform
y_tilde = y;
y_tilde(:,2:4) = log(y(:,2:4));

yTest = rand( size(y_tilde,1), size(y_tilde,2)-1);
thetaTest = rand(size(pram_true.theta));
testDynamicalModel(@hes1logmodelODE, @hes1logmodelDx, @hes1logmodelDtheta, ... 
                   "Hes1 log", yTest, thetaTest, y(:,1))
               
hes1model.fOde = @hes1logmodelODE;
hes1model.fOdeDx = @hes1logmodelDx;
hes1model.fOdeDtheta= @hes1logmodelDtheta;
hes1model.thetaLowerBound = [0 0 0 0 0 0 0];
hes1model.thetaUpperBound = [Inf Inf Inf Inf Inf Inf Inf];

config.sigma = pram_true.sigma;
config.useFixedSigma =  true;

hes1result = MagiSolver( y_tilde, hes1model, [], config);

% Trace plots of parameters and log-post
theta_names = ["a", "b", "c", "d", "e", "f", "g"];
plotMagiOutput(hes1result, "trace", theta_names, true, true, "mean", 0.025, 0.975, false, true, 4)

% Parameter estimates
summaryMagiOutput(hes1result, "mean", false, theta_names, 0.025, 0.975)

% Inferred trajectories (log scale)
comp_names = ["P (17 observations)", "M (16 observations)", "H (unobserved)"];
plotMagiOutput(hes1result, "traj", comp_names, true, true, "mean", 0.025, 0.975, [], [], 3)

% Inferred trajectories (original scale)
xEst = exp(squeeze(mean(hes1result.xsampled, 1)));
for i=1:size(xEst,2)
    subplot(1,3,i);
    qlim = zeros(2, size(hes1result.xsampled,2));
    for j=1:size(xEst,1)
        qlim(:,j) = exp(quantile(hes1result.xsampled(:,j,i),[0.025 0.975]));
    end
    fill([y(:,1)' fliplr(y(:,1)')],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(y(:,1)', y(:,i+1),'o', 'Color', 'black', 'linewidth', 2);
    plot(y(:,1)', xEst(:,i), 'LineWidth', 2, 'Color', 'green');
    plot(times, xtrue(:,i), 'LineWidth', 2, 'Color', 'red');
    xlabel('time');
    hold off;
    if i==3
        legend('95% cred. interval', 'noisy observations', 'inferred trajectory', 'truth',...
            'orientation', 'vertical', 'location', 'northeast');
    end
    title(comp_names(i));
end

save("hes1result.mat")


%%% Fitzhugh-Nagumo equations

clear;
FNdat = readmatrix('../data/FN-sim.csv');
rng(12321);

% Use setDiscretizationInterval to create evenly spaced set
y_I0 = setDiscretizationInterval(FNdat, 0.5);

% Use setDiscretization to insert additional evenly-spaced points between 
% each adjacent time point
y_I1 = setDiscretization(y_I0, 1);
y_I2 = setDiscretization(y_I0, 2);
y_I3 = setDiscretization(y_I0, 3);

fnmodel.fOde = @fnmodelODE;
fnmodel.fOdeDx = @fnmodelDx;
fnmodel.fOdeDtheta= @fnmodelDtheta;
fnmodel.thetaLowerBound= [0 0 0];
fnmodel.thetaUpperBound= [Inf Inf Inf];

config.niterHmc = 10000;

FNres0 = MagiSolver( y_I0, fnmodel, [], config);
FNres1 = MagiSolver( y_I1, fnmodel, [], config);
FNres2 = MagiSolver( y_I2, fnmodel, [], config);

config.nstepsHmc = 1000;
FNres3 = MagiSolver( y_I3, fnmodel, [], config);

% Visualize parameter estimates over different discretization sets
resList = {FNres0, FNres1, FNres2, FNres3};
theta_names = ["a", "b", "c", "sigmaV", "sigmaR"];
FNsummary = {};
for i=1:length(resList)
    FNsummary{i} = summaryMagiOutput(resList{i}, "mean", true, theta_names, 0.025, 0.975);
end

for i=1:length(theta_names)
    subplot(1,5,i);
    getMean = @(x) table2array(x(1,theta_names(i)));
    getErrlower = @(x) table2array(x(1,theta_names(i)))-table2array(x(2,theta_names(i)));
    getErrupper = @(x) table2array(x(3,theta_names(i)))-table2array(x(1,theta_names(i)));
    
    errorbar(1:4, cellfun(getMean, FNsummary), cellfun(getErrlower, FNsummary), cellfun(getErrupper, FNsummary),'o');
    xlim([0,5]);
    xticks(1:4);
    xticklabels({'I0', 'I1', 'I2', 'I3'});
    title(theta_names(i));
end

tvec = 0:0.01:20;

xEst_I0 = squeeze(mean(FNres0.xsampled, 1));
thetaEst_I0 = mean(FNres0.theta);

xEst_I1 = squeeze(mean(FNres1.xsampled, 1));
thetaEst_I1 = mean(FNres1.theta);

xEst_I2 = squeeze(mean(FNres2.xsampled, 1));
thetaEst_I2 = mean(FNres2.theta);

xEst_I3 = squeeze(mean(FNres3.xsampled, 1));
thetaEst_I3 = mean(FNres3.theta);

% Trajectories implied by estimates of parameters and initial conditions
[~, FNtr_I0]=ode45(@fnmodelODEsolve,tvec,xEst_I0(1,:),[],thetaEst_I0);
[~, FNtr_I1]=ode45(@fnmodelODEsolve,tvec,xEst_I1(1,:),[],thetaEst_I1);
[~, FNtr_I2]=ode45(@fnmodelODEsolve,tvec,xEst_I2(1,:),[],thetaEst_I2);
[~, FNtr_I3]=ode45(@fnmodelODEsolve,tvec,xEst_I3(1,:),[],thetaEst_I3);

subplot(1,2,1);
plot( tvec, horzcat(FNtr_I0(:,1), FNtr_I1(:,1), FNtr_I2(:,1), FNtr_I3(:,1)));
legend('I0','I1','I2','I3', 'orientation', 'horizontal', 'location', 'north', ...
    'AutoUpdate', 'off');
xlabel('Time');
title('V');
hold on;
plot( FNdat(:,1), FNdat(:,2), 'o', 'Color', 'black');
hold off;

subplot(1,2,2);
plot( tvec, horzcat(FNtr_I0(:,2), FNtr_I1(:,2), FNtr_I2(:,2), FNtr_I3(:,2)));
legend('I0','I1','I2','I3', 'orientation', 'horizontal', 'location', 'north', ...
    'AutoUpdate', 'off');
xlabel('Time');
title('R');
hold on;
plot( FNdat(:,1), FNdat(:,3), 'o', 'Color', 'black');
hold off;

% RMSDs
FNrmsd = vertcat(sqrt(mean((FNtr_I0(ismember(tvec, FNdat(:,1)),:) - FNdat(:,2:3)).^2,1)),...
    sqrt(mean((FNtr_I1(ismember(tvec, FNdat(:,1)),:) - FNdat(:,2:3)).^2,1)),...
    sqrt(mean((FNtr_I2(ismember(tvec, FNdat(:,1)),:) - FNdat(:,2:3)).^2,1)),...
    sqrt(mean((FNtr_I3(ismember(tvec, FNdat(:,1)),:) - FNdat(:,2:3)).^2,1)))';
array2table(FNrmsd,'VariableNames',{'I0','I1','I2','I3'}, 'RowNames', {'V','R'})

save("FNresult.mat")


%%% HIV time-dependent model

clear;
hivdtmodel.fOde = @hivtdmodelODE;
hivdtmodel.fOdeDx = @hivtdmodelDx;
hivdtmodel.fOdeDtheta= @hivtdmodelDtheta;
hivdtmodel.thetaLowerBound= [0 0 0 0 0];
hivdtmodel.thetaUpperBound= [Inf Inf Inf Inf Inf];

pram_true.theta = [36, 0.108, 0.5, 1000, 3]; % lambda, rho, delta, N, c
pram_true.x0 = [600, 30, 1e5]; % TU, TI, V
pram_true.sigma = [sqrt(10) sqrt(10) 10];
tvec = 0:0.2:20;

rng(12321);

[foo, xtrue]=ode45(@hivtdmodelODEsolve,tvec,pram_true.x0,[],pram_true.theta);
xtrue = horzcat(foo,xtrue); 

y = xtrue;
for j=1:(size(y,2)-1)
  y(:,1+j) = normrnd(y(:,1+j), pram_true.sigma(j));
end

compnames = ["TU", "TI", "V"];
complabels = ["Concentration", "Concentration", "Load"];

% use gpsmoothing to determine phi/sigma
phiEst = zeros(2, size(y,2)-1);
sigmaInit = zeros(1, size(y,2)-1);

for j=1:(size(y,2)-1)
    hyperparam = gpsmoothing( y(:,j+1), y(:,1), [], []);
    phiEst(:,j) = hyperparam.phi;
    sigmaInit(j) = hyperparam.sigma;
end
    
tOut = 0:0.05:20;
for i=1:3
  subplot(1,3,i);
  
  fitMean = gpmean(y(:,i+1), tvec, tOut, phiEst(:,i), sigmaInit(i), [], []);
  fitCov = gpcov(y(:,i+1), tvec, tOut, phiEst(:,i), sigmaInit(i), []);

  gp_UB = fitMean + 1.96 * sqrt(diag(fitCov));
  gp_LB = fitMean - 1.96 * sqrt(diag(fitCov));
  
  fill([tOut fliplr(tOut)],[gp_LB' fliplr(gp_UB')],[.9 .9 .9],'LineStyle','none')
  hold on;
  plot(tvec, y(:,i+1), '.', 'Color', 'black');

  xlabel('Time');
  ylabel(complabels(i));
  title(compnames(i));

  plot(tOut, fitMean, 'Color', "#0072BD");
  
  if (i<3)
      hold off;
  end
end

% override phi/sigma for V (3rd) component
phiEst(:,3) = [1e7, 0.5];
sigmaInit(3) = 100;

i = 3;
fitMean = gpmean(y(:,i+1), tvec, tOut, phiEst(:,i), sigmaInit(i), [], []);
fitCov = gpcov(y(:,i+1), tvec, tOut, phiEst(:,i), sigmaInit(i), []);

gp_UB = fitMean + 1.96 * sqrt(diag(fitCov));
gp_LB = fitMean - 1.96 * sqrt(diag(fitCov));

fill([tOut fliplr(tOut)],[gp_LB' fliplr(gp_UB')],[.9 .9 .9],'LineStyle','none');
plot(tvec, y(:,i+1), '.', 'Color', 'black');
plot(tOut, fitMean, 'Color', "#EDB120");
hold off;

config.sigma = sigmaInit;
config.phi = phiEst;

y_I = setDiscretization(y, 1);
HIVresult = MagiSolver( y_I, hivdtmodel, [], config);

% Parameter estimates
theta_names = ["lambda", "rho", "delta", "N", "c"];
summaryMagiOutput(HIVresult, "mean", false, theta_names, 0.025, 0.975)

% Inferred trajectories
xEst = squeeze(mean(HIVresult.xsampled, 1));
for i=1:size(xEst,2)
    subplot(1,3,i);
    qlim = zeros(2, size(HIVresult.xsampled,2));
    for j=1:size(xEst,1)
        qlim(:,j) = quantile(HIVresult.xsampled(:,j,i),[0.025 0.975]);
    end
    fill([y_I(:,1)' fliplr(y_I(:,1)')],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    %hold on;plot(y_I(:,1)', y_I(:,i+1),'.', 'Color', 'black', 'linewidth', 2);
    hold on;
    plot(y_I(:,1)', xEst(:,i), 'LineWidth', 2, 'Color', 'green');
    plot(tvec, xtrue(:,i+1), 'LineWidth', 1.5, 'Color', 'red');
    xlabel('Time');
    ylabel(complabels(i));
    hold off;
    if i==3
        legend('95% cred. interval', 'inferred trajectory', 'truth',...
            'orientation', 'vertical', 'location', 'northeast');
    end
    title(compnames(i));
end

save("HIVtdresult.mat")


%%% Hamiltonian Monte Carlo - example of sticky samples with high autocorrelation

load("FNresult.mat")

rng(12321);

y_I0 = setDiscretizationInterval(FNdat, 0.5);
y_I3 = setDiscretization(y_I0, 3);

config.niterHmc = 10000;
config.nstepsHmc = 200;
FNres3b = MagiSolver( y_I3, fnmodel, [], config);

theta_names = ["a", "b", "c", "sigmaV", "sigmaR"];
plotMagiOutput(FNres3b, "trace", theta_names, true, true, "mean", 0.025, 0.975, true, true, 3)

