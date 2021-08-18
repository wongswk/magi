function resultDx = HIVmodelDx(theta,x,t) 
  resultDx= zeros(size(x,1),size(x,2),size(x,2));

  T = exp(x(:,1));
  Tm = exp(x(:,2));
  Tw = exp(x(:,3));
  Tmw = exp(x(:,4));

  resultDx(:,1,1) = (0);
  resultDx(:,2,1) = -1e-6*theta(2)*Tm;
  resultDx(:,3,1) = -1e-6*theta(3)*Tw;
  resultDx(:,4,1) = -1e-6*theta(4)*Tmw;

  resultDx(:,1,2) = 1e-6*theta(2)*T + 1e-6*0.25*theta(4)*Tmw.*T./Tm;
  resultDx(:,2,2) = -1e-6*0.25*theta(4)*Tmw.*T ./ Tm;
  resultDx(:,3,2) = -1e-6*theta(5)*Tw;
  resultDx(:,4,2) = 0.25*1e-6*theta(4)*Tmw.*T./Tm;

  resultDx(:,1,3) = 1e-6*theta(3)*T + 0.25*1e-6*theta(4)*Tmw.*T./Tw;
  resultDx(:,2,3) = -1e-6*theta(6)*Tm;
  resultDx(:,3,3) = -1e-6*0.25*theta(4)*Tmw.*T ./ Tw;
  resultDx(:,4,3) = 1e-6*0.25*theta(4)*Tmw.*T./Tw;

  resultDx(:,1,4) = 1e-6*0.5*theta(4)*T;
  resultDx(:,2,4) = (1e-6*theta(5)+1e-6*theta(6))*Tw.*Tm./Tmw;
  resultDx(:,3,4) = (1e-6*theta(5)+1e-6*theta(6))*Tm.*Tw./Tmw;
  resultDx(:,4,4) = -(1e-6*theta(5)+1e-6*theta(6))*Tw.*Tm./Tmw;

end