function out = MagiSolver(y, odeModel, tvec, control)

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
    
    if isfield(control,'thetaInit') && ~isempty(control.thetaInit)
        thetaInitExogenous = control.thetaInit;
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
        stepSizeFactor = 0.01;
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


    [outCpp, phiUsed] = solveMagi(y, odeModel, tvec, sigmaExogenous, phiExogenous, xInitExogenous, thetaInitExogenous, ...
        muExogenous, dotmuExogenous, priorTemperatureLevel, priorTemperatureDeriv, 1, char('generalMatern'), nstepsHmc, ...
        burninRatio, niterHmc, stepSizeFactor, 1, bandSize, true, true, true, false, useFixedSigma, true);

    burnin = round(niterHmc*burninRatio);
    nI = size(y,1); % # points in I
    D = size(y,2); % # components

    out.lp = outCpp(1,(burnin+1):niterHmc)';
    out.xsampled = reshape(outCpp(2:(1+D*nI),(burnin+1):niterHmc)', [length((burnin+1):niterHmc) nI D]);
    out.theta = outCpp((2+D*nI):((2+(D)*nI)+length(odeModel.thetaLowerBound)-1),(burnin+1):niterHmc)';
    out.sigma = outCpp((2+D*nI+length(odeModel.thetaLowerBound)):size(outCpp,1),(burnin+1):niterHmc)';    
    out.phi = phiUsed;

end