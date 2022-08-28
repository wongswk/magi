function resultDtheta = HIVmodelDtheta(theta,x,tvec) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  T = exp(x(:,1));
  Tm = exp(x(:,2));
  Tw = exp(x(:,3));
  Tmw = exp(x(:,4));


  resultDtheta(:,1,1) = (1.0);
  resultDtheta(:,2,1) = -1e-6*Tm;

  resultDtheta(:,3,1) = -1e-6*Tw;

  resultDtheta(:,4,1) = -1e-6*Tmw;

  resultDtheta(:,2,2) = 1e-6*T;
  resultDtheta(:,4,2) = 1e-6*0.25*Tmw.*T ./ Tm;
  resultDtheta(:,5,2) = -1e-6*Tw;
  resultDtheta(:,7,2) = 1.0;

  resultDtheta(:,3,3) = 1e-6*T;
  resultDtheta(:,4,3) = 1e-6*0.25*Tmw.*T ./ Tw;
  resultDtheta(:,6,3) = -1e-6*Tm;
  resultDtheta(:,8,3) = 1.0;

  resultDtheta(:,4,4) = 1e-6*0.5 * T;
  resultDtheta(:,5,4) = 1e-6*Tw .* Tm ./ Tmw;
  resultDtheta(:,6,4) = 1e-6*Tw .* Tm ./ Tmw;
  resultDtheta(:,9,4) = 1.0;

end
