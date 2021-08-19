function resultDtheta = hes1logmodelDtheta(theta,x,tvec) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  P = x(:,1);
  M = x(:,2);
  H = x(:,3);


  resultDtheta(:,1,1) = -exp(H);
  resultDtheta(:,2,1) = exp(M-P);
  resultDtheta(:,3,1) = (-1);

  resultDtheta(:,4,2) = (-1);
  resultDtheta(:,5,2) = exp(-M)./(1+exp(2*P));

  resultDtheta(:,1,3) = -exp(P);
  resultDtheta(:,6,3) = exp(-H)./(1+exp(2*P));
  resultDtheta(:,7,3) = (-1);

end