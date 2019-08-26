function [value, gradient] = loglik(x)
    value = -0.5 * sum(x.^2);
    gradient = -x;
end
