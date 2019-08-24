function resultDtheta = fnmodelDtheta(theta,x) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  V = x(:,1);
  R = x(:,2);

  resultDtheta(:,3,1) = V - V.^3 / 3.0 + R;

  resultDtheta(:,1,2) =  1.0 / theta(3);

  resultDtheta(:,2,2) = -R / theta(3);
  resultDtheta(:,3,2) = 1.0/(theta(3)^2) * ( V - theta(1) + theta(2) * R);

end