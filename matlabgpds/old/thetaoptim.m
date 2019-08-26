function [f,  g] = thetaoptim(yobs, curCov, cursigma, theta, modelName, useBand )

  [f , g ] = xthetallikM( yobs, curCov, cursigma, theta, modelName, useBand );
  
  f = -f;
  g =- g  ( (size(yobs,1)*size(yobs,2)+1):length(theta));

end

