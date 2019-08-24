function resultDx = ptransmodelDx(theta,x) 
  resultDx = zeros(size(x,1),size(x,2),size(x,2));

  S = x(:,1);
  dS = x(:,2);
  R = x(:,3);
  RS = x(:,4);
  RPP = x(:,5);

  resultDx(:,1,1) = -theta(1) - theta(2) * R;
  resultDx(:,3,1) = -theta(2) * S;
  resultDx(:,4,1) = (theta(3));

  resultDx(:,1,2) = (theta(1));

  resultDx(:,1,3) = -theta(2)*R;
  resultDx(:,3,3) = -theta(2)*S;
  resultDx(:,4,3) = (theta(3));
  resultDx(:,5,3) =  theta(5) * theta(6) ./  (theta(6) + RPP).^2;

  resultDx(:,1,4) = theta(2)*R;
  resultDx(:,3,4) = theta(2)*S;
  resultDx(:,4,4) = (-theta(3) - theta(4));

  resultDx(:,4,5) = (theta(4));
  resultDx(:,5,5) = -theta(5) * theta(6) ./  (theta(6) + RPP).^2;

end
