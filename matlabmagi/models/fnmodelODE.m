function result = fnmodelODE(theta,x,tvec)
  V = x(:,1);
  R = x(:,2);
  
  result = zeros( size(x,1), size(x,2));
  result(:,1) = theta(3) * (V - V.^3 / 3.0 + R);
  result(:,2) = -1.0/theta(3) * ( V - theta(1) + theta(2) * R);
  
end