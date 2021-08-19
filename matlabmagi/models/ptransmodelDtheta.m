function resultDtheta = ptransmodelDtheta(theta,x,tvec) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  S = x(:,1);
  dS = x(:,2);
  R = x(:,3);
  RS = x(:,4);
  RPP = x(:,5);


  resultDtheta(:,1,1) = -S;
  resultDtheta(:,2,1) = -S.*R;
  resultDtheta(:,3,1) = RS;

  resultDtheta(:,1,2) = S;

  resultDtheta(:,2,3) = -S.*R;
  resultDtheta(:,3,3) = RS;
  resultDtheta(:,5,3) = RPP ./ (theta(6)+RPP);
  resultDtheta(:,6,3) = -theta(5) * RPP ./ (theta(6)+RPP).^2;

  resultDtheta(:,2,4) = S.*R;
  resultDtheta(:,3,4) = -RS;
  resultDtheta(:,4,4) = -RS;

  resultDtheta(:,4,5) = RS;
  resultDtheta(:,5,5) = - RPP ./ (theta(6)+RPP);
  resultDtheta(:,6,5) = theta(5) * RPP ./ (theta(6)+RPP).^2;
  
end