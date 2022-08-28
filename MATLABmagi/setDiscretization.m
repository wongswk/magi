function filled = setDiscretization(data, level)
%
% Set the discretization level of a data matrix for input to MagiSolver, by
% inserting time points where the GP is constrained to the derivatives of the ODE system.
% 
%      dat: data matrix. First column is time points.
%      level: discretization level (a positive integer). 2^level - 1 equally-spaced points will be inserted between existing data points in dat.
% 
% RETURN a data matrix with the same columns as dat, with rows added for the inserted discretization time points.

  if (level==0)
    filled = data;
    return
  end
  newdata = data;
  newdata = sortrows(newdata);
  dummydata = newdata(2:size(newdata,1),:);
  dummydata(:,:) = NaN;
  dummydata(:,1) = (newdata(2:size(newdata,1),1) + newdata(1:(size(newdata,1)-1),1))/2;
  newdata = vertcat(newdata, dummydata);
  newdata = sortrows(newdata);
  filled = setDiscretization(newdata, level-1);
