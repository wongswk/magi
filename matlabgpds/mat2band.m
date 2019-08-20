function [ A ] = mat2band( a, bandsize )
  N = size(a,1);
  A = zeros(2*bandsize+1,N);
  
  for j = 1:N
    k = bandsize + 1 - j;
    for i=max(1,j-bandsize):min(N,j+bandsize)
      A(k+i,j) = a(i,j);
    end
  end
  
  %A

end

