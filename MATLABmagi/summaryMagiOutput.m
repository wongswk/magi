function parEst = summaryMagiOutput(x, sigma, par_names, lower, upper)
%
% Computes a summary table of parameter estimates from the output of MagiSolver
% 
%      x: output of MagiSolver
%      sigma: logical; if true, the noise levels will be included in the summary.
%      par_names: vector of parameter names for the summary table. Should be the same length as the number of parameters in theta, or the combined length of theta and sigma when sigma is true.
%      lower: the lower quantile of the credible interval.
%      upper: the upper quantile of the credible interval.
%
% RETURN a table where rows display the posterior mean, lower credible limit, and upper credible limit of each parameter.
% 

    if (sigma)
        resMat = horzcat(x.theta, x.sigma);
    else
        resMat = x.theta;
    end
    parEst = vertcat(mean(resMat),quantile(resMat, [lower, upper]));
    parEst = array2table(parEst,'RowNames',{'Mean',[num2str(100*lower), '%'],[num2str(100*upper),'%']}, 'VariableNames', par_names);