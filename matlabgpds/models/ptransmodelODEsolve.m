function dx = ptransmodelODEsolve(t,x,pars)
    S = x(1);
    dS = x(2);
    R = x(3);
    RS = x(4);
    RPP = x(5);
    dx=x;

    dx(1) = -pars(1)*S - pars(2) * S .* R + pars(3) * RS;
    dx(2)= pars(1)*S;
    dx(3) = -pars(2)*S.*R + pars(3)*RS + pars(5) * RPP ./ (pars(6)+RPP);
    dx(4) = pars(2)*S.*R - pars(3)* RS - pars(4)*RS;
    dx(5) = pars(4)*RS - pars(5) * RPP ./ (pars(6)+RPP);
end