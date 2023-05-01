function out = MagiSolver(y, odeModel, tvec, control)
% MAnifold-constrained Gaussian process Inference (MAGI)
%
% Core function of the MAGI method for inferring the parameters and trajectories of 
% dynamic systems governed by ordinary differential equations.  See detailed usage examples included.
% 
%      y: data matrix of observations
%      odeModel: list of ODE functions and inputs. See details.
%      tvec: vector of discretization time points corresponding to rows of "y".
%            If missing, MagiSolver will use the first column of "y".
%      control: list of control variables, which may include `sigma`, `phi`, `xInit`,
%               `theta`, `mu`, `dotmu`, `priorTemperature`, `niterHmc`, `burninRatio`,
%               `nstepsHmc`, `stepSizeFactor`, `bandSize`, `useFixedSigma`,
%               `kerneltype`, `skipMissingComponentOptimization`, `positiveSystem`, `verbose`. See details.
% 
% RETURN
% MagiSolver returns a struct with the following elements:
%     theta: matrix of MCMC samples for the system parameters theta, after burn-in.
%     xsampled: array of MCMC samples for the system trajectories at each discretization time point, after burn-in.
%     sigma: matrix of MCMC samples for the observation noise SDs sigma, after burn-in.
%     phi: matrix of estimated GP hyper-parameters, one column for each system component.
%     lp: vector of log-posterior values at each MCMC iteration, after burn-in.
%     y, tvec, odeModel: the data matrix, time vector, and odeModel specification from the inputs to MagiSolver.
% 
% DETAILS
% The data matrix "y" has a column for each system component, and optionally with first column
% with the discretization time points. If the time column is not provided in `y`, a vector of time
% points must be provided via the `tvec` argument. The rows of `y` correspond to the discretization set
% I at which the GP is constrained to the derivatives of the ODE system. To set the desired discretization
% level for inference, use setDiscretization to prepare the data matrix for input into MagiSolver. 
% Missing observations are indicated with NaN.
% 
% The struct `odeModel` is used for specification of the ODE system and its parameters. It must include five elements:
%     fOde:  function that computes the ODEs, specified with the form f(theta, x, tvec).
%     fOdeDx:  function that computes the gradients of the ODEs with respect to the system components
%     fOdeDtheta:  function that computes the gradients of the ODEs with respect to the parameters theta
%     thetaLowerBound:  a vector indicating the lower bounds of each parameter in theta.
%     thetaUpperBound:  a vector indicating the upper bounds of each parameter in theta.
% 
% Additional control variables can be supplied to MagiSolver via the optional struct `control`, which may include the following:
%       sigma:  a vector of noise levels (observation noise standard deviations) sigma for each component, at which to initialize MCMC sampling.  
%               By default, MagiSolver computes starting values for sigma via Gaussian process (GP) smoothing. If the noise levels are known, specify sigma together with useFixedSigma = true.
%       phi:  a matrix of GP hyper-parameters for each component, with rows for the kernel hyper-parameters and columns for the system components.
%             By default, MagiSolver estimates phi via an optimization routine.
%       theta:  a vector of starting values for the parameters theta, at which to initialize MCMC sampling. 
%               By default, MagiSolver uses an optimization routine to obtain starting values.
%       xInit:  a matrix of values for the system trajectories of the same dimension as y, at which to initialize MCMC sampling. 
%               Default is linear interpolation between the observed (non-missing) values of y and an optimization routine for entirely unobserved components of y.
%       mu:  a matrix of values for the mean function of the GP prior, of the same dimension as y. Default is a zero mean function.
%       dotmu:  a matrix of values for the derivatives of the GP prior mean function, of the same dimension as y. Default is zero.
%       priorTemperature:  the tempering factor by which to divide the contribution of the GP prior, to control the influence of the GP prior relative to the likelihood.
%                          Default is the total number of observations divided by the total number of discretization points.
%       niterHmc:  MCMC sampling from the posterior is carried out via Hamiltonian Monte Carlo (HMC). niterHmc specifies the number of HMC iterations to run.  Default is 20000 HMC iterations.
%       nstepsHmc:  the number of leapfrog steps per HMC iteration. Default is 200.
%       burninRatio:  the proportion of HMC iterations to be discarded as burn-in. Default is 0.5, which discards the first half of the MCMC samples.
%       stepSizeFactor:  initial leapfrog step size factor for HMC. Can be a specified as a scalar (applied to all posterior dimensions) or a vector (with length corresponding to the dimension of the posterior). Default is 0.01, and the leapfrog step size is automatically tuned during burn-in to achieve an acceptance rate between 60-90\%.
%       bandSize:  a band matrix approximation is used to speed up matrix operations, with default band size 20. Can be increased if MagiSolver returns an error indicating numerical instability.
%       useFixedSigma:  logical, set to true if sigma is known.  If useFixedSigma is true, the known values of sigma must be supplied via the sigma control variable.
%       kerneltype: the GP covariance kernel, `generalMatern` is the default and recommended choice. Other available choices are `matern`, `rbf`, `compact1`, `periodicMatern`.
%       skipMissingComponentOptimization: logical, set to true to skip automatic optimization for missing components. If skipMissingComponentOptimization = true, values for xInit and phi must be supplied for all system components. Default is false.
%       positiveSystem: logical, set to true if the system cannot be negative. Default is false.       
%       verbose: logical, set to true to output diagnostic and progress messages to the console.  
% 
% 
% REFERENCE
% Shihao Yang, Samuel WK Wong, SC Kou (2021). Inference of dynamic systems from noisy and sparse data via manifold-constrained Gaussian processes.
%   Proceedings of the National Academy of Sciences, 118 (15), e2020397118.
%
    if isempty(tvec)
       tvec = y(:,1)';
       disp('tvec is not specified, first column in y will be used as time');
       y(:,1) = [];
    end
    
    if ~isrow(tvec)
        tvec = tvec';
    end

    if isfield(control,'sigma') && ~isempty(control.sigma)
        sigmaExogenous = control.sigma;
        sigmaExogenous(~isfinite(sigmaExogenous)) = 1;
    else
        sigmaExogenous = [];
    end

    if isfield(control,'phi') && ~isempty(control.phi)
        phiExogenous = control.phi;
    else
        phiExogenous = [];
    end
    
    if isfield(control,'xInit') && ~isempty(control.xInit)
        xInitExogenous = control.xInit;
    else
        xInitExogenous = [];
    end
    
    if isfield(control,'theta') && ~isempty(control.theta)
        thetaInitExogenous = control.theta;
    else
        thetaInitExogenous = [];
    end
    
    if isfield(control,'mu') && ~isempty(control.mu)
        muExogenous = control.mu;
    else
        muExogenous = [];
    end

    if isfield(control,'dotmu') && ~isempty(control.dotmu)
        dotmuExogenous = control.dotmu;
    else
        dotmuExogenous = [];
    end

    if isfield(control,'priorTemperature') && ~isempty(control.priorTemperature)
        priorTemperatureLevel = control.priorTemperature;
        priorTemperatureDeriv = control.priorTemperature;
    else 
        priorTemperatureLevel = 1/mean(isfinite(y), 'all');
        priorTemperatureDeriv = 1/mean(isfinite(y), 'all');
    end
    %disp(priorTemperatureLevel);

    if isfield(control,'niterHmc') && ~isempty(control.niterHmc)
        niterHmc = control.niterHmc;
    else
        niterHmc = 20000;
    end
    
    if isfield(control,'burninRatio') && ~isempty(control.burninRatio)
        burninRatio = control.burninRatio;
    else
        burninRatio = 0.5;
    end

    if isfield(control,'nstepsHmc') && ~isempty(control.nstepsHmc)
        nstepsHmc = control.nstepsHmc;
    else
        nstepsHmc = 200;
    end

    if isfield(control,'stepSizeFactor') && ~isempty(control.stepSizeFactor)
        stepSizeFactor = control.stepSizeFactor;
    else
        stepSizeFactor = [];
    end    

    if isfield(control,'bandSize') && ~isempty(control.bandSize)
        bandSize = control.bandSize;
    else
        bandSize = 20;
    end      
    
    if isfield(control,'useFixedSigma') && ~isempty(control.useFixedSigma)
        useFixedSigma = control.useFixedSigma;
    else
        useFixedSigma = false;
    end      

    if isfield(control,'kerneltype') && ~isempty(control.kerneltype)
        kerneltype = char(control.kerneltype);
    else
        kerneltype = char('generalMatern');
    end        
    
    if isfield(control,'skipMissingComponentOptimization') && ~isempty(control.skipMissingComponentOptimization)
        skipMissingComponentOptimization = control.skipMissingComponentOptimization;
    else
        skipMissingComponentOptimization = false;
    end    

    if isfield(control,'positiveSystem') && ~isempty(control.positiveSystem)
        positiveSystem = control.positiveSystem;
    else
        positiveSystem = false;
    end    
    
    if isfield(control,'verbose') && ~isempty(control.verbose)
        verbose = control.verbose;
    else
        verbose = false;
    end



    [outCpp, phiUsed] = solveMagi(y, odeModel, tvec, sigmaExogenous, phiExogenous, xInitExogenous, thetaInitExogenous, ...
        muExogenous, dotmuExogenous, priorTemperatureLevel, priorTemperatureDeriv, 1, kerneltype, nstepsHmc, ...
        burninRatio, niterHmc, stepSizeFactor, 1, bandSize, true, true, true, false, useFixedSigma, skipMissingComponentOptimization, positiveSystem, verbose);

    burnin = round(niterHmc*burninRatio);
    nI = size(y,1); % # points in I
    D = size(y,2); % # components

    out.lp = outCpp(1,(burnin+1):niterHmc)';
    out.xsampled = reshape(outCpp(2:(1+D*nI),(burnin+1):niterHmc)', [length((burnin+1):niterHmc) nI D]);
    out.theta = outCpp((2+D*nI):((2+(D)*nI)+length(odeModel.thetaLowerBound)-1),(burnin+1):niterHmc)';
    out.sigma = outCpp((2+D*nI+length(odeModel.thetaLowerBound)):size(outCpp,1),(burnin+1):niterHmc)';    
    out.phi = phiUsed;
    out.y = y;
    out.tvec = tvec;
    out.odeModel = odeModel;

end
