function plotMagiOutput(x, obs, ci, comp_names, lower, upper)
%
% Plots the inferred trajectories from the output of MagiSolver
% 
%      x: output of MagiSolver
%      obs: logical; if true, points will be added on the plots for the observations.
%      ci: logical; if true, credible bands will be added to the plots.
%      comp_names: vector of system component names. Should be the same
%      length as the number of system components.
%      lower: the lower quantile of the credible band.
%      upper: the upper quantile of the credible band.
%
% RETURN plots of inferred trajectories (posterior means) and credible bands from the MCMC samples, one subplot for each system component. 
% 

    xEst = squeeze(mean(x.xsampled, 1));
    for i=1:size(xEst,2)
        subplot(1,size(xEst,2),i);
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
        plot(x.tvec', xEst(:,i), 'LineWidth', 2, 'Color', 'green');
        xlabel('Time');
        hold off;
        title(comp_names(i));
    end
