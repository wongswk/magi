function [ ret ] = calCov( phi, rInput, signrInput, bandsize, noiseInjection, complexity)
  if nargin < 4 || isempty(bandsize)
    bandsize = size(rInput,1);
  end
  
  if nargin < 5 || isempty(noiseInjection)
      noiseInjection = 1e-7;
  end
  
  if nargin < 6
      complexity = 3;
  end
  
  

  ret = calCovGeneralMatern(phi, rInput, signrInput);
  
  d1 = size(ret.C,1);
  
  ret.mu = zeros(d1);
  ret.dotmu = zeros(d1);
  ret.C = ret.C + noiseInjection * eye( size(rInput,1));
  
  if complexity < 2
      return
  end
    
  ret.dCdphiCube = zeros(d1,d1,2);
  ret.dCdphiCube(:,:,1) = ret.Cdphi1;
  ret.dCdphiCube(:,:,2) = ret.Cdphi2;
  
  [ret.CeigenVec, ret.Ceigen1over] = eig(ret.C);
  ret.Ceigen1over = 1./diag(ret.Ceigen1over);
  ret.Cinv = ret.CeigenVec* ( ret.Ceigen1over .* transpose(ret.CeigenVec) );
  ret.mphi = ret.Cprime * ret.Cinv;
  
  ret.Kright = sqrt(ret.Ceigen1over) .* transpose(ret.CeigenVec) * transpose(ret.Cprime);
  ret.Kphi = ret.Cdoubleprime - transpose(ret.Kright) * ret.Kright  + noiseInjection .* eye( size(rInput,1));
  
  [ret.KeigenVec, ret.Keigen1over] = eig(ret.Kphi);
  ret.Keigen1over = 1./diag(ret.Keigen1over);
  ret.Kinv = ret.KeigenVec * (ret.Keigen1over.*transpose(ret.KeigenVec));
  
  
%     Keigen1over <- 1/Kdecomp$values
%     KeigenVec <- Kdecomp$vectors
%     Kinv <- KeigenVec%*%(Keigen1over*t(KeigenVec))
%         
%     mphiLeftHalf <- Cprime %*% CeigenVec  
  
%   retmore <- with(ret, {
%     dCdphiCube <- sapply(dCdphi, identity, simplify = "array")
%     Cdecomp <- eigen(C)
%     Ceigen1over <- 1/Cdecomp$value
%     CeigenVec <- Cdecomp$vectors
%     Cinv <- CeigenVec%*%(Ceigen1over*t(CeigenVec))
%     mphi <-  Cprime %*% Cinv
%     
%     Kright <- sqrt(Ceigen1over) * t(CeigenVec) %*% t(Cprime)
%     Kphi <- Cdoubleprime - t(Kright)%*%Kright  + noiseInjection * diag( nrow(rInput))
%     Kdecomp <- eigen(Kphi)
%     Keigen1over <- 1/Kdecomp$values
%     KeigenVec <- Kdecomp$vectors
%     Kinv <- KeigenVec%*%(Keigen1over*t(KeigenVec))
%         
%     mphiLeftHalf <- Cprime %*% CeigenVec
%     list(Ceigen1over = Ceigen1over,
%          CeigenVec = CeigenVec,
%          Cinv = Cinv, 
%          mphi = mphi, 
%          Kphi = Kphi, 
%          Keigen1over = Keigen1over,
%          KeigenVec = KeigenVec,
%          Kinv = Kinv,
%          mphiLeftHalf = mphiLeftHalf,
%          dCdphiCube = dCdphiCube)
%   })
  %retmore <- bandCov(retmore, bandsize)
  ret.CinvBand = mat2band(ret.Cinv, bandsize);
  ret.mphiBand = mat2band(ret.mphi, bandsize);
  ret.KinvBand = mat2band(ret.Kinv, bandsize);
  ret.bandsize = bandsize;
  
  ret.Cdphi1 = [];
  ret.Cdphi2 = [];
  
%   out <- c(ret, retmore)
%   out$dCdphi <- NULL
%   out

end

