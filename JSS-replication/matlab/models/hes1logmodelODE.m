function PMHdt = hes1logmodelODE(theta,x,tvec) 
  P = exp(x(:,1));
  M = exp(x(:,2));
  H = exp(x(:,3));


  PMHdt= zeros(size(x,1),size(x,2));
  PMHdt(:,1) = -theta(1)*H + theta(2)*M./P - theta(3);
  PMHdt(:,2) = -theta(4) + theta(5)./(1+P.^2)./M;
  PMHdt(:,3) = -theta(1)*P + theta(6)./(1+P.^2)./H - theta(7);

end