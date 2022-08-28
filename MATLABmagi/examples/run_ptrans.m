addpath('models/');

% Ensure libcmagi.so or libcmagi.dll is accessible to MATLAB

config.nobs = 15;
config.t_end = 100;
config.linfillspace = 0.5;
config.bandSize = 40;
config.ndis = config.t_end / config.linfillspace + 1;
config.priorTemperature = config.ndis / config.nobs;

pram_true.theta = [0.07,0.6,0.05,0.3,0.017,0.3];
pram_true.x0 = [1,0,1,0,0];

% observation time points
tobs = [0,1,2,4,5,7,10,15,20,30,40,50,60,80,100]; 
% example noisy observation data simulated from this system
dat = [0.999079656108803,0.586308085200401,0.416974630597961,0.247679615712861,0.184309297707158,0.119944831062102,0.0681865620595117,0.0236552286911344,0.00745687186276209,0.00992390748438715,-0.00576258115637723,0.0113716708361056,0.00712444177328236,0.00338374623734491,-0.00647203743140396;
-0.0206688222850907,0.0743013677255439,0.0875508322569361,0.131818149826651,0.143629966415659,0.161945662971894,0.192283137514025,0.201993834320625,0.194551102073354,0.193769661527693,0.224513380314124,0.210511172809517,0.192373842438558,0.193916763266096,0.207562851574146;
1.00189800441379,0.629359915424173,0.498304359515507,0.379228373564295,0.344017232712093,0.34053743984727,0.345075538313128,0.349626747470537,0.405795575839972,0.528840197511724,0.628097608053192,0.704614931258064,0.794672928576815,0.902737985033323,0.962546553407464;
0.00756816885217754,0.307035997193477,0.354267874665295,0.28415022038115,0.243440228129352,0.16647935033977,0.0744357980673708,0.0344091660323895,0.0169200158004037,0.0149706273296933,-0.00584039336162013,-0.0273533741101581,-0.00212029807624876,0.0181129218246689,0.00084728445535421;
0.00797818050672776,0.0638545616266391,0.138781291478125,0.357998833530703,0.379614221210997,0.504191787764218,0.584246911770709,0.605418521607663,0.583254256902994,0.485336221460864,0.389136506748854,0.297568678380025,0.22450821683515,0.102600678389572,0.0448000401086414];

xsim_obs = horzcat(tobs', dat');

% To simulate from system using another random seed, can do the following:
% config.noise = 0.01 * ones(1,5); % 0.01 high noise, 0.001 low noise
% config.seed = rand(1)*1e7;
% times = 0:0.1:config.t_end;
% [foo, xtrue]=ode45(@ptransmodelODEsolve,times,pram_true.x0,[],pram_true.theta);
% xtrue = horzcat(foo,xtrue); 
% 
% xsim = tobs;
% for j=2:size(xtrue,2)
%     xsim = horzcat(xsim,interp1(times,xtrue(:,j),xsim(:,1)));
% end
% 
% rng(config.seed);
% for j=2:size(xsim,2)
%   xsim(:,j) = normrnd(xsim(:,j), config.noise(j-1));
% end
% 
% xsim_obs = xsim;
%%%%%%%%

xsim = setDiscretizationInterval(xsim_obs, config.linfillspace);

xsimt = xsim(:,1);
for j=2:size(xsim,2) % linear interpolate starting values
    xsimt = horzcat(xsimt,interp1(xsim_obs(:,1),xsim_obs(:,j),xsimt(:,1)));
end
xInit = xsimt(:,2:size(xsim,2));

% inference ----------------------------
ptmodel.fOde = @ptransmodelODE;
ptmodel.fOdeDx = @ptransmodelDx;
ptmodel.fOdeDtheta= @ptransmodelDtheta;
ptmodel.thetaLowerBound= [0 0 0 0 0 0];
ptmodel.thetaUpperBound= [4 4 4 4 4 4];

% get initial phi and sigma estimate at set I0
config.niterHmc = 2;
res1 = MagiSolver( xInit(ismember(xsim(:,1),0:1:100),:), ptmodel, 0:1:100, config);

% re-optimize phi
config.sigma = res1.sigma;
res2 = MagiSolver( xInit(ismember(xsim(:,1),0:1:100),:), ptmodel, 0:1:100, config);

% HMC sampling
config.sigma = res2.sigma;
config.phi = res2.phi;
config.xInit = xInit;
config.niterHmc = 20000;
gpode = MagiSolver( xsim, ptmodel, [], config);

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
[foo, xRecons]=ode45(@ptransmodelODEsolve,times,xEst(1,:),[],thetaEst);
for i=1:size(xEst,2)
    subplot(1,size(xEst,2),i);
    plot(times,xRecons(:,i));
    hold on;
    plot(xsim_obs(:,1)', xsim_obs(:,i+1),'Marker','*','LineStyle','none');
    hold off;
end
