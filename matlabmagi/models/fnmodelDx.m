function resultDx = fnmodelDx(theta,x,t) 
  resultDx = zeros(size(x,1),size(x,2),size(x,2));

  V = x(:,1);

  resultDx(:,1,1) = theta(3) * (1 - V.^2);
  resultDx(:,2,1) = theta(3);

  resultDx(:,1,2) = (-1.0 / theta(3));
  resultDx(:,2,2) = ( -1.0*theta(2)/theta(3) );

end