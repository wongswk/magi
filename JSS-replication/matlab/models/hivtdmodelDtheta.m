function resultDtheta = hivtdmodelDtheta(theta,x,tvec) 
  TU = x(:,1);
  TI = x(:,2);
  V = x(:,3);

  lambda = theta(1);
  rho = theta(2);
  delta = theta(3);
  N = theta(4);
  c = theta(5);

  eta = 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000));
  
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));
  
  resultDtheta(:,1,1) = 1;
  resultDtheta(:,2,1) = -TU;
  resultDtheta(:,3,1) = 0;
  resultDtheta(:,4,1) = 0;
  resultDtheta(:,5,1) = 0;

  resultDtheta(:,1,2) = 0;
  resultDtheta(:,2,2) = 0;
  resultDtheta(:,3,2) = -TI;
  resultDtheta(:,4,2) = 0;
  resultDtheta(:,5,2) = 0;

  resultDtheta(:,1,3) = 0;
  resultDtheta(:,2,3) = 0;
  resultDtheta(:,3,3) = N * TI;
  resultDtheta(:,4,3) = delta * TI;
  resultDtheta(:,5,3) = -V;