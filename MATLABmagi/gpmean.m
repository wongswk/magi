function ret = gpmean(yobs, tvec, tnew, phi, sigma, kerneltype, deriv)
% Conditional mean of Gaussian process given observations
%
% Compute the conditional mean of a Gaussian process (and optionally, its derivative), given a vector of observations, hyper-parameters phi, and noise standard deviation sigma.
% 
%      yobs: vector of observations
%      tvec: vector of time points corresponding to observations
%      tnew: vector of time points at which the conditional mean should be computed
%      phi: vector of hyper-parameters for the covariance kernel (kerneltype)
%      sigma: the noise standard deviation.
%      kerneltype: the covariance kernel, types `matern`, `compact1`,
%           `periodicMatern`, `generalMatern` are supported. Default (if blank) is `generalMatern`.
%      deriv: logical; if true, the conditional mean of the GP's derivative is also computed
% 
% RETURN A vector with the values of the conditional mean function evaluated at the time points in tnew.
% If deriv = true, returned as a struct where .mean is the mean and .deriv contains the values of the conditional mean of the GP derivative evaluated at the time points in tnew.

    if isempty(kerneltype)
        kerneltype = char('generalMatern');
    else
        kerneltype = char(kerneltype);
    end
    
    if isempty(deriv)
        deriv = false;
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

    
    res = calcMeanCurve(tvec, yobs, tnew, phi, sigma, kerneltype, deriv);
    
    if (deriv)
        ret.mean = res(:,:,1);
        ret.deriv = res(:,:,2);
    else        
        ret = res(:,:,1);
    end


