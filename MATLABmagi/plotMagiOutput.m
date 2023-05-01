function plotMagiOutput(x, type, comp_names, obs, ci, est, lower, upper, sigma, lp, nplotcol)
%
% Plots the inferred trajectories or diagnostic traceplots from the output of MagiSolver
% 
%      x: output of MagiSolver
%      type: string; type = "traj" plots inferred trajectories, while setting type = "trace" generates diagnostic traceplots for the MCMC samples of the parameters and log-posterior values.
%      comp_names: vector of names for the plot panels. Should be the same length as the number of system components (when type = "traj") or parameters (when type = "trace").
%      obs: logical; if true, points will be added on the plots for the observations when type = "traj".
%      ci: logical; if true, credible bands/intervals will be added to the plots.
%      est: string specifying the posterior quantity to plot as the estimate. Can be "mean", "median", "mode", or "none". Recommended is "mean", which plots the posterior mean of the MCMC samples.
%      lower: the lower quantile of the credible band.
%      upper: the upper quantile of the credible band.
%      sigma: logical; if true, the noise levels \sigmaσ will be included in the traceplots when type = "trace".
%      lp: logical; if true, the values of the log-posterior will be included in the traceplots when type = "trace".
%      nplotcol: the number of subplots per row.
%
% RETURN Plots the inferred system trajectories (when type = "traj") or diagnostic traceplots of the parameters and log-posterior (when type = "trace") from the MCMC samples. The posterior mean is recommended as the estimate of the trajectories and parameters (est = "mean"). Alternatives are the posterior median (est = "median", taken component-wise) and the posterior mode (est = "mode", approximated by the MCMC sample with the highest log-posterior value).
% Setting type = "traj" produces plots of the inferred trajectories and credible bands from the MCMC samples, one subplot for each system component; lower = 0.025 and upper = 0.975 produces a central 95% credible band when ci = TRUE. Adding the observed data points (obs = TRUE) can provide a visual assessment of the inferred trajectories.
% Setting type = "trace" generates diagnostic traceplots for the MCMC samples of the system parameters and the values of the log-posterior, which is a useful tool for informally assessing convergence. In this case, the est and ci options add horizontal lines to the plots that indicate the estimate (in red) and credible interval (in green) for each parameter.
% 
    if (est == "mode")
        lpmaxInd = find(x.lp == max(x.lp));
    end

    if (type == "traj")
        
        xEst = squeeze(mean(x.xsampled, 1));
        nplotrow = ceil(size(xEst, 2)/nplotcol);
        if (size(xEst,2) < nplotcol)
            nplotcol = size(xEst,2);
        end
        
        for i=1:size(xEst,2)
            subplot(nplotrow,nplotcol,i);
            
            qlim = zeros(2, size(x.xsampled,2));
            for j=1:size(xEst,1)
                qlim(:,j) = quantile(x.xsampled(:,j,i),[lower upper]);
            end
            if (ci)
                fill([x.tvec fliplr(x.tvec)],[(qlim(1,:)) (fliplr(qlim(2,:)))],[.9 .9 .9],'LineStyle','none')
            end
            hold on;
            if (obs)
                plot(x.tvec', x.y(:,i), 'o', 'Color', 'black', 'linewidth', 2);
            end
            
            if (est == "mean")
                plot(x.tvec', xEst(:,i), 'LineWidth', 2, 'Color', 'green');
            end
            if (est == "median")
                xMed = squeeze(median(x.xsampled, 1));
                plot(x.tvec', xMed(:,i), 'LineWidth', 2, 'Color', 'green');
            end
            if (est == "mode")
                plot(x.tvec', x.xsampled(lpmaxInd,:,i), 'LineWidth', 2, 'Color', 'green');
            end
            
            xlabel('Time');
            hold off;
            title(comp_names(i));
        end
    end
    
    if (type == "trace")
        
        if (lp)
            comp_names = [comp_names, "log-post"];
        end
        
        allpar = x.theta;
        if (sigma)
            allpar = horzcat(allpar, x.sigma);
        end
        if (lp)
            allpar = horzcat(allpar, x.lp);
        end
        
        nplotrow = ceil(size(allpar, 2)/nplotcol);
        if (size(allpar,2) < nplotcol)
            nplotcol = size(allpar,2);
        end        
        
        for i=1:size(allpar,2)
            subplot(nplotrow,nplotcol,i);
            plot(allpar(:,i), 'Color', 'black');
            hold on;
            
            if (est == "mean")
                yline( mean(allpar(:,i)), 'Color', 'red');
            end
            if (est == "median")
                yline( median(allpar(:,i)), 'Color', 'red');
            end
            if (est == "mode")
                yline( allpar(lpmaxInd,i), 'Color', 'red');
            end
            
            if (ci)
                yline(quantile(allpar(:,i),lower), 'Color', 'green');
                yline(quantile(allpar(:,i),upper), 'Color', 'green');
            end
            
            hold off;
            title(comp_names(i));
        end
        
        
    end
            
        
