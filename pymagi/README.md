# Usage in Python

### Installation

Compile the Python MAGI library by running `py_build.sh`. Before running the examples, ensure that this `pymagi` directory is within Python's path.

### Examples

Four complete, self-contained dynamic system examples are provided in the `examples/` directory.

  * Hes1 oscillation system
  * FitzHugh-Nagumo (FN) system
  * Protein transduction (PTrans) system
  * HIV time-dependent system

Each example lays out how to prepare the input to `MagiSolver`, and then use the output to plot inferred trajectories and compute parameter estimates.

Further details on how to set up your own ODE system for inference are provided below, following the FitzHugh-Nagumo system  as an illustrative example.

### Preparing the system of differential equations

First, we create a function that defines the ODEs.  It takes as input a parameter vector `theta`, matrix `x` (where each column of `x` is a system component), and time vector `tvec`: 

```
def fOde(theta, x, tvec):
    V = x[:, 0]
    R = x[:, 1]
    Vdt = theta[2] * (V - pow(V,3) / 3.0 + R)
    Rdt = -1.0/theta[2] * ( V - theta[0] + theta[1] * R)
    result = np.stack([Vdt, Rdt], axis=1)
    return result
```

Next, we create functions that calculate the gradient of the ODEs, (1) with respect to the components `x`, and (2) with respect to the parameters `theta`, that provide output as 3-D arrays as follows:

```
def fOdeDx(theta, x, tvec):
    resultDx = np.zeros(shape=[np.shape(x)[0], np.shape(x)[1], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDx[:,0,0] = theta[2] * (1 - np.square(V))
    resultDx[:,1,0] = theta[2]
    resultDx[:,0,1] = -1.0 / theta[2]
    resultDx[:,1,1] = -1.0*theta[1]/theta[2]
    return resultDx

def fOdeDtheta(theta, x, tvec):
    resultDtheta = np.zeros(shape=[np.shape(x)[0], np.shape(theta)[0], np.shape(x)[1]])
    V = x[:, 0]
    R = x[:, 1]
    resultDtheta[:,2,0] = V - pow(V,3) / 3.0 + R
    resultDtheta[:,0,1] = 1.0 / theta[2]
    resultDtheta[:,1,1] = -R / theta[2]
    resultDtheta[:,2,1] = 1.0/pow(theta[2], 2) * ( V - theta[0] + theta[1] * R)
    return resultDtheta
```
Now we are ready to create an `ode_system` object representing the dynamic system model to pass to `MagiSolver`.  First, give the system a name.  Then, supply the three arguments -- ODEs, gradient with respect to `x`, gradient with respect to `theta` -- using the three functions created above.  Lastly, supply two vectors, `thetaLowerBound` and `thetaUpperBound`, that specify the lower and upper bounds on the parameters `theta`.

```
fn_system = ode_system("FN-python", fOde, fOdeDx, fOdeDtheta,
                       thetaLowerBound=np.array([0,0,0]), thetaUpperBound=np.array([np.inf, np.inf, np.inf]))

```








