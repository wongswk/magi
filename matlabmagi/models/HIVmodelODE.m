function HIVdt = HIVmodelODE(theta,x,t) 
  T = exp(x(:,1));
  Tm = exp(x(:,2));
  Tw = exp(x(:,3));
  Tmw = exp(x(:,4));

  HIVdt= zeros(size(x,1),size(x,2));

  HIVdt(:,1) = (theta(1) - 1e-6*theta(2)*Tm - 1e-6*theta(3)*Tw - 1e-6*theta(4)*Tmw);
  HIVdt(:,2) = (theta(7) + 1e-6*theta(2)*T - 1e-6*theta(5)*Tw) + 1e-6*0.25*theta(4)*Tmw.*T ./ Tm;
  HIVdt(:,3) = (theta(8) + 1e-6*theta(3)*T - 1e-6*theta(6)*Tm) + 1e-6*0.25*theta(4)*Tmw.*T ./ Tw;
  HIVdt(:,4) = theta(9) + 0.5*1e-6*theta(4)*T + (1e-6*theta(5)+1e-6*theta(6))*Tw.*Tm ./ Tmw;

end
