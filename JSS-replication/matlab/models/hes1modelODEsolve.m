function dx = hes1modelODEsolve(t,x,theta) 
  P = x(1);
  M = x(2);
  H = x(3);
  dx = x;

  dx(1) = -theta(1)*P*H + theta(2)*M - theta(3)*P;
  dx(2) = -theta(4)*M + theta(5)/(1+P^2);
  dx(3) = -theta(1)*P*H + theta(6)/(1+P^2) - theta(7)*H;

end
