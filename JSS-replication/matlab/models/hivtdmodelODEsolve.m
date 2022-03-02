function dx = hivtdmodelODEsolve(t,x,pars) 
  TU = x(1);
  TI = x(2);
  V = x(3);

  lambda = pars(1);
  rho = pars(2);
  delta = pars(3);
  N = pars(4);
  c = pars(5);

  eta = 9e-5 * (1 - 0.9 * cos(pi * t / 1000));

  dx = x;
  dx(1) = lambda - rho * TU - eta .* TU .* V;
  dx(2) = eta .* TU .* V - delta * TI;
  dx(3) = N * delta * TI - c * V;
