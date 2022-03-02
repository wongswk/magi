% MATLAB replication script

% To run this script, ensure that path contains 
% - libcmagi.so (or libcmagi.dll on Windows) 
% - 'mex' compiled solveMagi.cpp and gpsmooth.cpp
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
y(:,2:4) = log(y(:,2:4));

yTest = rand( size(y,1), size(y,2)-1);
thetaTest = rand(size(pram_true.theta));
testDynamicalModel(@hes1logmodelODE, @hes1logmodelDx, @hes1logmodelDtheta, ... 
                   "Hes1 log", yTest, thetaTest, y(:,1));
               
hes1model.fOde = @hes1logmodelODE;
hes1model.fOdeDx = @hes1logmodelDx;
hes1model.fOdeDtheta= @hes1logmodelDtheta;
hes1model.thetaLowerBound = [0 0 0 0 0 0 0];
hes1model.thetaUpperBound = [Inf Inf Inf Inf Inf Inf Inf];

config.sigma = pram_true.sigma;
config.useFixedSigma =  true;

hes1result = MagiSolver( y, hes1model, [], config);

% Trace plots of parameters and log-post
theta_names = ["a", "b", "c", "d", "e", "f", "g"];
for i=1:size(hes1result.theta,2)
    subplot(2, 4, i);
    plot(hes1result.theta(:,i));
    title(theta_names(i));
end
subplot(2, 4, 8);
plot(hes1result.lp);
title("log-post");

% Parameter estimates
hes1est = vertcat(mean(hes1result.theta),quantile(hes1result.theta, [0.025, 0.975]));
array2table(hes1est,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

% Inferred trajectories
comp_names = ["P (17 observations)", "M (16 observations)", "H (unobserved)"];
xEst = exp(squeeze(mean(hes1result.xsampled, 1)));
for i=1:size(xEst,2)
    subplot(1,3,i);
    qlim = zeros(2, size(hes1result.xsampled,2));
    for j=1:size(xEst,1)
        qlim(:,j) = exp(quantile(hes1result.xsampled(:,j,i),[0.025 0.975]));
    end
    fill([y(:,1)' fliplr(y(:,1)')],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(y(:,1)', exp(y(:,i+1)),'o', 'Color', 'black', 'linewidth', 2);
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

% Parameter estimates
theta_names = ["a", "b", "c", "sigmaV", "sigmaR"];

FNest0 = vertcat(mean(horzcat(FNres0.theta, FNres0.sigma)), ...
    quantile(horzcat(FNres0.theta, FNres0.sigma), [0.025, 0.975]));
array2table(FNest0,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

FNest1 = vertcat(mean(horzcat(FNres1.theta, FNres1.sigma)), ...
    quantile(horzcat(FNres1.theta, FNres1.sigma), [0.025, 0.975]));
array2table(FNest1,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

FNest2 = vertcat(mean(horzcat(FNres2.theta, FNres2.sigma)), ...
    quantile(horzcat(FNres2.theta, FNres2.sigma), [0.025, 0.975]));
array2table(FNest2,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

FNest3 = vertcat(mean(horzcat(FNres3.theta, FNres3.sigma)), ...
    quantile(horzcat(FNres3.theta, FNres3.sigma), [0.025, 0.975]));
array2table(FNest3,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

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

for i=1:3
  subplot(1,3,i);
  plot(tvec, xtrue(:,i+1));
  xlabel('Time');
  ylabel(complabels(i));
  title(compnames(i));
  hold on;
  plot(tvec, y(:,i+1), '.');
end

% use gpsmoothing to determine phi/sigma
phiEst = zeros(2, size(y,2)-1);
sigmaInit = zeros(1, size(y,2)-1);

for j=1:(size(y,2)-1)
    hyperparam = gpsmoothing( y(:,j+1), y(:,1), [], []);
    phiEst(:,j) = hyperparam.phi;
    sigmaInit(j) = hyperparam.sigma;
end
    
% override phi/sigma for V (3rd) component
phiEst(:,3) = [1e7, 0.5];
sigmaInit(3) = 100;

config.sigma = sigmaInit;
config.phi = phiEst;

y_I = setDiscretization(y, 1);
HIVresult = MagiSolver( y_I, hivdtmodel, [], config);

% Parameter estimates
theta_names = ["lambda", "rho", "delta", "N", "c"];
HIVest = vertcat(mean(HIVresult.theta),quantile(HIVresult.theta, [0.025, 0.975]));
array2table(HIVest,'RowNames',{'Mean','2.5%','97.5%'}, 'VariableNames', theta_names)

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
