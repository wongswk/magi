function PMHdt = hes1modelODE(theta,x,tvec) 
  P = x(:,1);
  M = x(:,2);
  H = x(:,3);


  PMHdt = zeros(size(x,1),size(x,2));
  PMHdt(:,1) = -theta(1)*P.*H + theta(2)*M - theta(3)*P;
  PMHdt(:,2) = -theta(4)*M + theta(5)./(1+P.^2);
  PMHdt(:,3) = -theta(1)*P.*H + theta(6)./(1+P.^2) - theta(7)*H;

end
