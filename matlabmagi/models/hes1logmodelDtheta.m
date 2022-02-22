function resultDtheta = hes1logmodelDtheta(theta,x,tvec) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  logP = x(:,1);
  logM = x(:,2);
  logH = x(:,3);


  resultDtheta(:,1,1) = -exp(logH);
  resultDtheta(:,2,1) = exp(logM-logP);
  resultDtheta(:,3,1) = (-1);

  resultDtheta(:,4,2) = (-1);
  resultDtheta(:,5,2) = exp(-logM)./(1+exp(2*logP));

  resultDtheta(:,1,3) = -exp(logP);
  resultDtheta(:,6,3) = exp(-logH)./(1+exp(2*logP));
  resultDtheta(:,7,3) = (-1);

end