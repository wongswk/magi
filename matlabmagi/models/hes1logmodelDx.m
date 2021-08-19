function resultDx = hes1logmodelDx(theta,x,tvec) 
  resultDx= zeros(size(x,1),size(x,2),size(x,2));

  P = x(:,1);
  M = x(:,2);
  H = x(:,3);


  expMminusP = exp(M-P);
  dP = -(1+exp(2*P)).^(-2).*exp(2*P)*2;

  resultDx(:,1,1) = -theta(2).*expMminusP;
  resultDx(:,2,1) = theta(2).*expMminusP;
  resultDx(:,3,1) = -theta(1)*exp(H);

  resultDx(:,1,2) = theta(5)*exp(-M).*dP;
  resultDx(:,2,2) = -theta(5)*exp(-M)./(1+exp(2*P));

  resultDx(:,1,3) = -theta(1)*exp(P) + theta(6)*exp(-H).*dP;
  resultDx(:,3,3) = -theta(6)*exp(-H)./(1+exp(2*P));

end