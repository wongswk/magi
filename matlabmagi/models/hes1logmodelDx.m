function resultDx = hes1logmodelDx(theta,x,tvec) 
  resultDx= zeros(size(x,1),size(x,2),size(x,2));

  logP = x(:,1);
  logM = x(:,2);
  logH = x(:,3);


  expMminusP = exp(logM-logP);
  dP = -(1+exp(2*logP)).^(-2).*exp(2*logP)*2;

  resultDx(:,1,1) = -theta(2).*expMminusP;
  resultDx(:,2,1) = theta(2).*expMminusP;
  resultDx(:,3,1) = -theta(1)*exp(logH);

  resultDx(:,1,2) = theta(5)*exp(-logM).*dP;
  resultDx(:,2,2) = -theta(5)*exp(-logM)./(1+exp(2*logP));

  resultDx(:,1,3) = -theta(1)*exp(logP) + theta(6)*exp(-logH).*dP;
  resultDx(:,3,3) = -theta(6)*exp(-logH)./(1+exp(2*logP));

end