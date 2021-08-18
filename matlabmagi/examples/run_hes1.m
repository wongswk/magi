addpath('models/');

% Ensure libcmagi.so or libcmagi.dll is accessible to MATLAB

config.nobs = 33;
config.noise = [0.15 0.15 NaN];
config.seed = rand(1)*1e7;
config.useFixedSigma =  true;
config.t_end=240;

pram_true.theta = [0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3];
pram_true.x0 = [log(1.438575), log(2.037488), log(17.90385)];
pram_true.sigma=config.noise;

times = 0:0.1:config.t_end;
[foo, xtrue]=ode45(@hes1logmodelODEsolve,times,pram_true.x0,[],pram_true.theta);

xtrue = horzcat(foo,xtrue); 

% Simulate noisy observations 
rng(config.seed);
xsim = xtrue;
for j=1:(size(xsim,2)-1)
  xsim(:,1+j) = normrnd(xsim(:,1+j), config.noise(j));
end

xsim_obs = xsim( linspace(1, size(xsim,1), config.nobs),:);

% Set asynchronous observation schedule
xsim_obs(:,4) = NaN;
xsim_obs(2:2:config.nobs,2) = NaN;
xsim_obs(1:2:config.nobs,3) = NaN;

xsim = setDiscretization(xsim_obs,0);

% inference ----------------------------
hes1model.fOde = @hes1logmodelODE;
hes1model.fOdeDx = @hes1logmodelDx;
hes1model.fOdeDtheta= @hes1logmodelDtheta;
hes1model.thetaLowerBound= [0 0 0 0 0 0 0];
hes1model.thetaUpperBound= [Inf Inf Inf Inf Inf Inf Inf];

gpode = MagiSolver( xsim, hes1model, [], config);

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
[foo, xRecons]=ode45(@hes1logmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:size(xEst,2)
    subplot(1,size(xEst,2),i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end