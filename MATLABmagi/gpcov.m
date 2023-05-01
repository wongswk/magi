function ret = gpcov(yobs, tvec, tnew, phi, sigma, kerneltype)
% Conditional covariance of Gaussian process given observations
%
% Compute the conditional covariance of a Gaussian process, given a vector of observations, hyper-parameters phi, and noise standard deviation sigma.
% 
%      yobs: vector of observations
%      tvec: vector of time points corresponding to observations
%      tnew: vector of time points at which the conditional covariance should be computed
%      phi: vector of hyper-parameters for the covariance kernel (kerneltype)
%      sigma: the noise standard deviation.
%      kerneltype: the covariance kernel, types `matern`, `compact1`,
%           `periodicMatern`, `generalMatern` are supported. Default (if blank) is `generalMatern`.
% 
% RETURN The conditional covariance matrix for the GP evaluated at the time points in tnew.

    if isempty(kerneltype)
        kerneltype = char('generalMatern');
    else
        kerneltype = char(kerneltype);
    end
    
    if ~isrow(yobs)
        yobs = yobs';
    end
    
    if ~isrow(tvec)
        tvec = tvec';
    end

    if ~isrow(tnew)
        tnew = tnew';
    end

    
    res = calcCovCurve(tvec, yobs, tnew, phi, sigma, kerneltype);
    
    ret = res(:,:,1);


