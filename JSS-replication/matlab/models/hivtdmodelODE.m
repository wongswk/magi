function result = hivtdmodelODE(theta,x,tvec) 
  TU = x(:,1);
  TI = x(:,2);
  V = x(:,3);

  lambda = theta(1);
  rho = theta(2);
  delta = theta(3);
  N = theta(4);
  c = theta(5);

  eta = 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000));

  result = zeros(size(x,1),size(x,2));
  result(:,1) = lambda - rho * TU - eta .* TU .* V;
  result(:,2) = eta .* TU .* V - delta * TI;
  result(:,3) = N * delta * TI - c * V;
