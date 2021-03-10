% Metropolis-Hastings using numerical solver

fn.ode=@fhn;
fn.like=@fhnlik;
fn.prior=@fhnprior;

reps = 50;
reppars = repmat(0,reps,4);
time=0:2:20; % observation times
n.iter=10000;
init=[-1 1];  % initial conditions V and R
pars=[.2,.2,3,.2]; % a,b,c,sigma

for j=1:reps
    [foo path]=ode45(fn.ode,time,init,[],pars(1:3));
    
    y = normrnd(path,pars(4));
    
    pars_start=[.2 .2 3 .2];  % start MCMC at truth
    mypars=repmat(0,n.iter,4);
    mypars(1,:)=pars_start;
    
    stepvar=[.1, .1, .1, .1]*0.1;
    accepts=0;
    
    tic;
    for k=2:(n.iter) % run M-H
        
        proposal= normrnd(mypars(k-1,:),stepvar);
        X=proposal;
        [loglike,path2]=fn.like(y,time,X,fn);
        log_alpha_numer = loglike+fn.prior(X,path2);
        
        [loglike,path]=fn.like(y,time,mypars(k-1,:),fn);
        log_alpha_denom = loglike+fn.prior(mypars(k-1,:),path);
        
        % make a decision
        if(rand<exp(log_alpha_numer - log_alpha_denom))
            accepts=accepts+1;
            mypars(k,:)=proposal;
        else
            mypars(k,:)=mypars(k-1,:);
            
            
        end
    end
    toc;
    
    %figure; for lp=1:4;subplot(2,4,lp);plot(mypars(1:n.iter,lp));end
    
    reppars(j,1:4) = mean(mypars((0.5*n.iter+1):n.iter,:)) % posterior mean after burn-in
    

end


