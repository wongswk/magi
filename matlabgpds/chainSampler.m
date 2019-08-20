function [ out ] = chainSampler(config, xInit, singleSampler, stepLowInit, verbose)

%chainSampler <- function(config, xInit, singleSampler, stepLowInit, verbose=TRUE, 
%                         thetaDim=3){
  numparam = length(xInit);
  n_iter = config.n_iter;
  xth_formal = zeros(n_iter, numparam);
  stepLow = stepLowInit;
  accepts = zeros(1, n_iter);
  accepts(1) = 0;
  lliklist =  zeros(1, n_iter);
  xth_formal(1,:) = xInit;
  burnin = round(n_iter*config.burninRatio);
  
  for t=2:n_iter
    rstep = stepLow + stepLow .* rand(1, length(stepLow));
    foo = singleSampler(xth_formal(t-1,:), rstep);
    
    xth_formal(t,:) = foo.final;
    accepts(t) = foo.acc;
    
    if t < burnin && t > 10
      if mean(accepts(max(1,t-99):t)) > 0.9
        stepLow = stepLow * 1.005;
      elseif mean(accepts(max(1,t-99):t)) < 0.6
        stepLow = stepLow * .995;
      end
      if mod(t, 100) == 0 
        xthsd = std(xth_formal(max(1,t-99):t,:));
        if mean(xthsd)>0
            stepLow = 0.05 .*xthsd/mean(xthsd)*mean(stepLow) + 0.95*stepLow;
        end
      end
    end
    lliklist(t) = foo.lpr;
    
    if verbose && mod(t, 100) == 0 
      fprintf(' Iter: %.0f, %.4f, %.4f \n', t, mean(accepts((t-99):t)), foo.lpr);
    end
    
  end
  %list(xth=xth.formal, lliklist=lliklist, stepLow=stepLow)

  out.xth = xth_formal;
  out.lliklist = lliklist;
  out.stepLow = stepLow;
end

