function parEst = summaryMagiOutput(x, est, sigma, par_names, lower, upper)
%
% Computes a summary table of parameter estimates from the output of MagiSolver
% 
%      x: output of MagiSolver
%      est: string specifying the posterior quantity to treat as the estimate. Recommend est = "mean", which treats the posterior mean as the estimate. Alternatives are the posterior median (est = "median", taken component-wise) and the posterior mode (est = "mode", approximated by the MCMC sample with the highest log-posterior value).
%      sigma: logical; if true, the noise levels will be included in the summary.
%      par_names: vector of parameter names for the summary table. Should be the same length as the number of parameters in theta, or the combined length of theta and sigma when sigma is true.
%      lower: the lower quantile of the credible interval.
%      upper: the upper quantile of the credible interval.
%
% RETURN a table where rows display the estimate, lower credible limit, and upper credible limit of each parameter.
% 

    
    if (est == "mean")
        f = @(x) mean(x);
        est_lab = "Mean";
    end
    if (est == "median")
        f = @(x) median(x);
        est_lab = "Median";
    end
    if (est == "mode")
        lpmaxInd = find(x.lp == max(x.lp));
        f = @(x) x(lpmaxInd,:);
        est_lab = "Mode";
    end
        

    if (sigma)
        resMat = horzcat(x.theta, x.sigma);
    else
        resMat = x.theta;
    end
    parEst = vertcat(f(resMat),quantile(resMat, [lower, upper]));
    parEst = array2table(parEst,'RowNames',{char(est_lab),[num2str(100*lower), '%'],[num2str(100*upper),'%']}, 'VariableNames', par_names);