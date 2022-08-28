function filled = setDiscretizationInterval(dat, by)
%
% Set the discretization level of a data matrix for input to MagiSolver, by
% inserting time points where the GP is constrained to the derivatives of the ODE system.
% 
%      dat: data matrix. First column is time points.
%      by: discretization interval. Equally-spaced spaced time points will be inserted with interval "by" between successive points.
% 
% RETURN a data matrix with the same columns as dat, with rows added for the inserted discretization time points.
    t_end = max(dat(:,1));
    fillC = 0:by:t_end;
    filled = zeros(length(fillC),size(dat,2));
    filled(:,:) = NaN;
    filled(:,1) = fillC';
    for i=1:size(filled,1)
        [isin, loc] = ismember(filled(i,1), dat(:,1));
        if isin
            filled(i,2:size(dat,2)) = dat(loc,2:size(dat,2));
        end
    end
