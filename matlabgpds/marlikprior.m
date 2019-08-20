function [f,g] = marlikprior(par, xsim_obs, r_nobs, config, j, priormeanFactor, priorsdFactor)
    [f, g] = phisigllikC( par, xsim_obs(:,1+j), r_nobs, char(config.kernel));
    f_penalty = log(normpdf(par(2), max(xsim_obs(:,1))*priormeanFactor,max(xsim_obs(:,1))*priorsdFactor));
    g_penalty = (par(2) - max(xsim_obs(:,1))*priormeanFactor) / (max(xsim_obs(:,1))*priorsdFactor)^2;
    f =  - (f + f_penalty);
    g = -g;
    g(2) = g(2) + g_penalty;
end