function [ logprior ] = fhnprior ( pars,path)
    a_val = log(normpdf( pars(1),0,.4));
    b_val = log(normpdf( pars(2),0,.4));
    c_val = log(chi2pdf( pars(3),2));
    sigma2_val = 1 / pars(4);
    %init_val = normpdf(init(1),0,2) * normpdf(init(2),0,2);
    
    logprior = a_val + b_val + c_val +  sigma2_val;

end

