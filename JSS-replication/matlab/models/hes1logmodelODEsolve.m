function dx = hes1logmodelODEsolve(t,x,theta) 
  P = exp(x(1));
  M = exp(x(2));
  H = exp(x(3));
  dx = x;

  dx(1) = -theta(1)*H + theta(2)*M./P - theta(3);
  dx(2) = -theta(4) + theta(5)./(1+P.^2)./M;
  dx(3) = -theta(1)*P + theta(6)./(1+P.^2)./H - theta(7);

end