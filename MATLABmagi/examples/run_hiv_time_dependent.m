addpath('models/');

% Ensure libcmagi.so or libcmagi.dll is accessible to MATLAB

config.nobs = 101;
config.noise = [sqrt(10) sqrt(10) 10];
config.seed = rand(1)*1e7;
config.t_end = 20;
config.filllevel = 0;

pram_true.theta = [36, 0.108, 0.5, 1000, 3]; % lambda, rho, delta, N, c
pram_true.x0 = [600, 30, 1e5]; % TU, TI, V
pram_true.sigma = config.noise;

times = 0:0.1:config.t_end;
[foo, xtrue]=ode45(@hivtdmodelODEsolve,times,pram_true.x0,[],pram_true.theta);
xtrue = horzcat(foo,xtrue); 

xsim  = linspace(0,config.t_end,config.nobs)';
for j=2:size(xtrue,2)
    xsim = horzcat(xsim,interp1(times,xtrue(:,j),xsim(:,1)));
end

rng(config.seed);
for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

xsim_obs = xsim;

xsim = setDiscretization(xsim_obs,config.filllevel);

testDynamicalModel(@hivtdmodelODE, @hivtdmodelDx, @hivtdmodelDtheta, "HIV time-dependent", xsim_obs(:,2:4), pram_true.theta, xsim_obs(:,1))

% inference ----------------------------
hivdtmodel.fOde = @hivtdmodelODE;
hivdtmodel.fOdeDx = @hivtdmodelDx;
hivdtmodel.fOdeDtheta= @hivtdmodelDtheta;
hivdtmodel.thetaLowerBound= [0 0 0 0 0];
hivdtmodel.thetaUpperBound= [Inf Inf Inf Inf Inf];

% use gpsmoothing to determine phi/sigma
phiExogenous = zeros(2, size(xsim,2)-1);
sigmaInit = zeros(1,size(xsim,2)-1);

for j=1:(size(xsim,2)-1)
    hyperparam = gpsmoothing( xsim_obs(:,j+1), xsim_obs(:,1), [], []);
    phiExogenous(:,j) = hyperparam.phi;
    sigmaInit(j) = hyperparam.sigma;
end
    
% override phi/sigma for V (3rd) component
phiExogenous(:,3) = [5e7, 1];
sigmaInit(3) = 1;

config.sigma = sigmaInit;
config.phi = phiExogenous;

gpode = MagiSolver( xsim, hivdtmodel, [], config);

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
[foo, xRecons]=ode45(@hivtdmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:size(xEst,2)
    subplot(1,size(xEst,2),i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end
