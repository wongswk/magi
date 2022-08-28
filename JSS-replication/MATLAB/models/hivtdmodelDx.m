function resultDx = hivtdmodelDx(theta,x,tvec) 
  TU = x(:,1);
  TI = x(:,2);
  V = x(:,3);

  lambda = theta(1);
  rho = theta(2);
  delta = theta(3);
  N = theta(4);
  c = theta(5);

  eta = 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000));

  resultDx= zeros(size(x,1),size(x,2),size(x,2));
  
  resultDx(:,1,1) = -rho - eta .* V;
  resultDx(:,2,1) = 0;
  resultDx(:,3,1) = -eta .* TU;

  resultDx(:,1,2) = eta .* V;
  resultDx(:,2,2) = -delta;
  resultDx(:,3,2) = eta .* TU;

  resultDx(:,1,3) = 0;
  resultDx(:,2,3) = N * delta;
  resultDx(:,3,3) = -c;