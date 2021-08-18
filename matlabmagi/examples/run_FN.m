addpath('models/');

% Ensure libcmagi.so or libcmagi.dll is accessible to MATLAB

config.nobs = 41;
config.noise = [0.2, 0.2];
config.filllevel = 2;
config.t_end = 20;

pram_true.theta = [0.2,0.2,3];
pram_true.x0 = [-1, 1];
pram_true.sigma=config.noise;

% observation time points
tobs = 0:0.5:20;
% sample noisy data
V = [ -1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, ...
  0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, ...
  -0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, ... 
  -1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84];
R = [0.94, 1.22, 0.89, 0.13, 0.4, 0.04, -0.21, -0.65, -0.31, -0.65,... 
  -0.72, -1.26, -0.56, -0.44, -0.63, 0.21, 1.07, 0.57, 0.85, 1.04, ...
  0.92, 0.47, 0.27, 0.16, -0.41, -0.6, -0.58, -0.54, -0.59, -1.15, ...
  -1.23, -0.37, -0.06, 0.16, 0.43, 0.73, 0.7, 1.37, 1.1, 0.85, 0.23];

xsim_obs = horzcat(tobs', V', R');

% To simulate from system using another random seed, can do the following:
% times = 0:0.1:config.t_end;
% [foo, xtrue]=ode45(@fnmodelODEsolve,times,pram_true.x0,[],pram_true.theta);
% 
% xtrue = horzcat(foo,xtrue); 
% 
% xsim  = linspace(0,config.t_end,config.nobs)';
% for j=2:size(xtrue,2)
%     xsim = horzcat(xsim,interp1(times,xtrue(:,j),xsim(:,1)));
% end
% 
% config.seed = rand(1)*1e7;
% rng(config.seed);
% 
% for j=1:(size(xsim,2)-1)
%   xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
% end
% 
% xsim_obs = xsim;
%%%%%%%

xsim = setDiscretization(xsim_obs,config.filllevel);

% inference ----------------------------
fnmodel.fOde = @fnmodelODE;
fnmodel.fOdeDx = @fnmodelDx;
fnmodel.fOdeDtheta= @fnmodelDtheta;
fnmodel.thetaLowerBound= [0 0 0];
fnmodel.thetaUpperBound= [Inf Inf Inf];

gpode = MagiSolver( xsim(:,2:size(xsim,2)), fnmodel, xsim(:,1)', config);

% Inferred trajectory and parameter estimates
xEst = squeeze(mean(gpode.xsampled, 1));
thetaEst = mean(gpode.theta);

% Sampled trajectories
for i=1:size(xEst,2)
    qlim = zeros(2, size(gpode.xsampled,2));
    for j=1:size(xEst,1)
        qlim(:,j) = quantile(gpode.xsampled(:,j,i),[0.025 0.975]);
    end

    subplot(1,size(xEst,2),i)
    fill([xsim(:,1)' fliplr(xsim(:,1)')],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
    hold on;plot(xsim(:,1)', (xsim(:,i+1)),'Marker','*')
    plot(xsim(:,1)', xEst(:,i), 'LineWidth', 1, 'Color', 'blue');
    hold off;
end

% Histograms of parameters
for i=1:size(gpode.theta,2)
    subplot(1, size(gpode.theta,2), i);
    histogram(gpode.theta(:,i));
    xline(pram_true.theta(i), 'LineWidth', 2, 'Color', 'red');
end
    
% Look at whether these estimates are reasonable for reconstructing trajectories using ODE solver
times = 0:0.1:config.t_end;
[foo, xRecons]=ode45(@fnmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:size(xEst,2)
    subplot(1,size(xEst,2),i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end

