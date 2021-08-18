function resultdt = ptransmodelODE(theta,x,t) 
  S = x(:,1);
  dS = x(:,2);
  R = x(:,3);
  RS = x(:,4);
  RPP = x(:,5);

  resultdt= zeros(size(x,1),size(x,2));

  resultdt(:,1) = -theta(1)*S - theta(2) * S .* R + theta(3) * RS;
  resultdt(:,2) = theta(1)*S;
  resultdt(:,3) = -theta(2)*S.*R + theta(3)*RS + theta(5) * RPP ./ (theta(6)+RPP);
  resultdt(:,4) = theta(2)*S.*R - theta(3)* RS - theta(4)*RS;
  resultdt(:,5) = theta(4)*RS - theta(5) * RPP ./ (theta(6)+RPP);

end
