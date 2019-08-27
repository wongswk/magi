function [ ydy ] = getMeanCurve( x, y, x_new, phi_mat, sigma_mat)
  %tvec <- c(x.new,x)
  
  %foo <- outer(tvec, t(tvec),'-')[,1,]
  %r <- abs(foo)
  %r2 <- r^2
  
  tvec = vertcat(x_new, x);
  
  foo1 = @(x,y)x-y;
  foo = bsxfun(foo1, tvec, tvec');
  r = abs(foo);
  r2 = r.^2;
  signr = -sign(foo);
  
  y_new = zeros(size(phi_mat,1), length(x_new));
  dy_new = zeros(size(phi_mat,1), length(x_new));
  
  for it=1:size(phi_mat,1)
    sigma = sigma_mat(it);
    phi = phi_mat(it,:)';
    
    covObj = calCov(phi, r, signr, [],[], 1);
    C = covObj.C;

    C((1+(size(C,1)+1)*length(x_new)):(size(C,1)+1):end) = C((1+(size(C,1)+1)*length(x_new)):(size(C,1)+1):end) + sigma^2;
    
    %diag(C)[-(1:length(x.new))] <- diag(C)[-(1:length(x.new))]+sigma^2
    y_new(it,:) = C(1:length(x_new), (length(x_new)+1):size(C,2)) * (C((length(x_new)+1):size(C,2), (length(x_new)+1):size(C,2)) \ y);    %* solve(C[-(1:length(x.new)),-(1:length(x.new))], y)
    %if(deriv){
    dy_new(it,:) = covObj.Cprime(1:length(x_new), 1:length(x_new)) * ( covObj.C(1:length(x_new), 1:length(x_new)) \ y_new(it,:)');
    %}
  end
  
  %if(deriv){
  %  return(list(y.new, dy.new))
  %}else{
  %  return(y.new)  
  %}
  ydy.y = y_new;
  ydy.dy = dy_new;
  

end

