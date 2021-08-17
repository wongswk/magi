function resultDtheta = hes1modelDtheta(theta,x,t) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  P = x(:,1);
  M = x(:,2);
  H = x(:,3);

  resultDtheta(:,1,1) = -P .* H;
  resultDtheta(:,2,1) = M;
  resultDtheta(:,3,1) = -P;

  resultDtheta(:,4,2) = -M;
  resultDtheta(:,5,2) = 1./(1 + P.^2);

  resultDtheta(:,1,3) = -P .* H;
  resultDtheta(:,6,3) = 1./(1 + P.^2);
  resultDtheta(:,7,3) = -H;

end
