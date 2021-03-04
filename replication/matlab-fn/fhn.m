function dx = fhn(t,x,pars) 
    a = pars(1);
    b = pars(2);
    c = pars(3);

	V=x(1);
	R=x(2);
    
    dx=x;

	dx(1) = c*(V-V^3/3+R);
	dx(2)= -1/c*(V - a+b*R);

	%dx = [dV;dR];
