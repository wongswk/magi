function [ covRet ] = calCovGeneralMatern( phi, r, signr, df)

  if nargin < 4
      df = 2.01;
  end
  %r2 = r.^2;
  
  x4bessel = sqrt(2.0 * df) .* r / phi(2);
  
  C = phi(1) * 2^(1-df) * exp(-gammaln(df)) * x4bessel.^df .* besselk(df,x4bessel);
  C(r==0) = phi(1);
  
  dCdphi2 = C .* (df ./ x4bessel - (besselk(df-1,x4bessel) + besselk(df+1,x4bessel))./(2*besselk(df,x4bessel)));
  dCdphi2 = dCdphi2 .* (-sqrt(2.0 * df) * r ./ phi(2)^2);
  dCdphi2(isnan(dCdphi2)) = 0;
  
  
  Cprime = C .* (df ./ x4bessel - (besselk(df-1,x4bessel) + besselk(df+1,x4bessel))./(2*besselk(df,x4bessel)));
  Cprime = Cprime .* sqrt(2.0 * df) .* signr / phi(2);
  Cprime(isnan(Cprime)) = 0;
  Cprime = -Cprime;
  
  Cdoubleprime = -phi(1) * 2^(1-df) * exp(-gammaln(df)) * 2.0 * df / phi(2)^2 * (df*(df-1).*x4bessel.^(df-2).*besselk(df,x4bessel) - df.*x4bessel.^(df-1).*(besselk(df-1,x4bessel)+besselk(df+1,x4bessel))...
    + x4bessel.^df.*(besselk(df-2,x4bessel)+2*besselk(df,x4bessel)+besselk(df+2,x4bessel))/4);

  
  Cdoubleprime(1:(size(Cdoubleprime,1)+1):end) = phi(1) * 2^(1-df) * exp(-gammaln(df)) * 2.0 * df / phi(2)^2 * gamma(df-1) * 2^(df-2);
  
  %dCdphi <- list(
  %  C/phi[1],
  %  dCdphi2
  %)
  %return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
  covRet.C = C;
  covRet.Cprime = Cprime;
  covRet.Cdoubleprime = Cdoubleprime;
  covRet.Cdphi1 = C/phi(1);
  covRet.Cdphi2 = dCdphi2;

end

