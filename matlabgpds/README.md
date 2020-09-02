# Usage in MATLAB

Three complete, self-contained dynamic system examples are provided in the `examples/` directory.

  * Hes1 oscillation system
  * FitzHugh-Nagumo (FN) system
  * Protein transduction (PTrans) system

Each example lays out how to prepare the input to `GpdsSolver`, and then use the output to plot inferred trajectories and compute parameter estimates.

Further details on how to set up your own ODE system for inference are provided below, following the FitzHugh-Nagumo system  as an illustrative example.

### Preparing the system of differential equations

First, we create a function that defines the ODEs.  It takes as input a parameter vector `theta` and matrix `x` where each column of `x` is a system component: 

```
function result = fnmodelODE(theta,x)
  V = x(:,1);
  R = x(:,2);
  
  result = zeros( size(x,1), size(x,2));
  result(:,1) = theta(3) * (V - V.^3 / 3.0 + R);
  result(:,2) = -1.0/theta(3) * ( V - theta(1) + theta(2) * R);
  
end
```

Next, we create functions that calculate the gradient of the ODEs, (1) with respect to the components `x`, and (2) with respect to the parameters `theta`, that provide output as 3-D arrays as follows:

```
function resultDx = fnmodelDx(theta,x) 
  resultDx = zeros(size(x,1),size(x,2),size(x,2));

  V = x(:,1);

  resultDx(:,1,1) = theta(3) * (1 - V.^2);
  resultDx(:,2,1) = theta(3);

  resultDx(:,1,2) = (-1.0 / theta(3));
  resultDx(:,2,2) = ( -1.0*theta(2)/theta(3) );

end

function resultDtheta = fnmodelDtheta(theta,x) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  V = x(:,1);
  R = x(:,2);

  resultDtheta(:,3,1) = V - V.^3 / 3.0 + R;

  resultDtheta(:,1,2) =  1.0 / theta(3);

  resultDtheta(:,2,2) = -R / theta(3);
  resultDtheta(:,3,2) = 1.0/(theta(3)^2) * ( V - theta(1) + theta(2) * R);

end
```
Now we are ready to create a `struct` representing the dynamic system model to pass to `GpdsSolver`.  Supply the three arguments -- ODEs, gradient with respect to `x`, gradient with respect to `theta` -- using the three functions created above.  Lastly, supply two vectors, `thetaLowerBound` and `thetaUpperBound`, that specify the lower and upper bounds on the parameters `theta`.

```
fnmodel.fOde = @fnmodelODE;
fnmodel.fOdeDx = @fnmodelDx;
fnmodel.fOdeDtheta= @fnmodelDtheta;
fnmodel.thetaLowerBound= [0 0 0];
fnmodel.thetaUpperBound= [Inf Inf Inf];
```








