function [ loglike,path] = fhnlik( data, time,pars,fn )

    sigma = pars(4);
    [time path] = ode15s(fn.ode,time,[-1 1],[],pars(1:3));
 %   d= size(data);
    % fprintf(1,'%d %d \n',d(1),d(2));
   %    fprintf(1,'%d %d',size(sol));
    if (sigma<=0)
        loglike = -1e20;
        return
    end
        
   
    %lik1 = sum( log(normpdf(data,sol,repmat(sqrt(sigma2),d(1),d(2)))));
    %like = lik1(1)+lik1(2);        
    %path=path(2:end,:);
    loglike=0;       
    for i=1:length(data)
      loglike=loglike+log(mvnpdf(data(i,:),path(i,:),sigma^2*eye(2)));
    end
        
    
end

