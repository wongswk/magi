function ret = gpsmoothing(yobs, tvec, kerneltype, sigma)
% Gaussian process smoothing
%
% Estimate hyper-parameters phi and noise standard deviation sigma for a vector of observations using Gaussian process smoothing.
% 
%      yobs: vector of observations
%      tvec: vector of time points corresponding to observations
%      kerneltype: the covariance kernel, types `matern`, `compact1`,
%           `periodicMatern`, `generalMatern` are supported. Default (if blank) is `generalMatern`.
%      sigma: the noise level (if known). By default, both phi and sigma are estimated.
%             If a value for sigma is supplied, then sigma is held fixed at the supplied value and only phi is estimated.
% 
% RETURN A struct containing phi and sigma with their estimated values.

    if isempty(kerneltype)
        kerneltype = char('generalMatern');
    else
        kerneltype = char(kerneltype);
    end
    
    if ~isrow(yobs)
        yobs = yobs';
    end    
    
    foo1 = @(x,y)x-y;
    foo = bsxfun(foo1, tvec, tvec');
    distInput = abs(foo);    
    yInput = yobs - mean(yobs);
    
    if (isempty(sigma))
        res = gpsmooth(yInput, distInput, kerneltype, -1, true);
        ret.sigma = res(end);
        ret.phi = res(1:(length(res)-1));
    else
        res = gpsmooth(yInput, distInput, kerneltype, sigma, true);
        ret.sigma = sigma;
        ret.phi = res;
    end
            


