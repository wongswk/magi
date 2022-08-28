function passed = testDynamicalModel(modelODE, modelDx, modelDtheta, modelName, x, theta, tvec)
% Given functions for the ODE and its gradients (with respect to the system components and parameters), 
% verify the correctness of the gradients using numerical differentiation.
%
%      modelODE: function that computes the ODEs
%      modelDx: function that computes the gradients of the ODEs with respect to the system components
%      modelDtheta: function that computes the gradients of the ODEs with respect to the parameters theta
%      modelName: string giving a name for the model
%      x: data matrix of system values, one column for each component, at which to test the gradients
%      theta: vector of parameter values for theta, at which to test the gradients
%      tvec: vector of time points corresponding to the rows of x
% 
% RETURN A struct with fields Dx and Dtheta, with value 1 if the corresponding gradient check passed and 0 if not.
%  
    deltaSmall = 1e-06;
    tolerance = 1e-04;
    passed.modelName = modelName;
    
    f = modelODE(theta, x, tvec);
    fdX = modelDx(theta, x, tvec);
    
    numericalDx = NaN( length(tvec), size(x,2), size(x,2));
    for j=1:size(x,2)
        xnew = x;
        xnew(:,j) = xnew(:,j) + deltaSmall;
        numericalDx(:,j,:) = (modelODE(theta, xnew, tvec) - f) / deltaSmall;
    end
    
    passed.Dx = max( abs(numericalDx - fdX), [], 'all') < tolerance;
    
    fDtheta = modelDtheta(theta, x, tvec);
    
    numericalDtheta = NaN( length(tvec), length(theta), size(x,2));
    for j=1:length(theta)
        thetanew = theta;
        thetanew(:,j) = thetanew(:,j) + deltaSmall;
        numericalDtheta(:,j,:) = (modelODE(thetanew, x, tvec) - f) / deltaSmall;
    end
    
    passed.Dtheta = max( abs(numericalDtheta - fDtheta), [], 'all') < tolerance;
    
    
    
    
    
    
        
    
    