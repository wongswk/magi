function resultDx = hes1modelDx(theta,x,t) 
  resultDx= zeros(size(x,1),size(x,2),size(x,2));

  P = x(:,1);
  H = x(:,3);

  resultDx(:,1,1) = -theta(1)*H - theta(3);
  resultDx(:,2,1) = ( theta(2) );
  resultDx(:,3,1) = -theta(1)*P;

  resultDx(:,1,2) = -2*theta(5)*P ./ (1.0 + P.^2).^2;
  resultDx(:,2,2) = ( -theta(4) );

  resultDx(:,1,3) = -theta(1)*H - 2*theta(6)*P ./ (1.0 + P.^2).^2;
  resultDx(:,3,3) = -theta(1)*P - theta(7);

end