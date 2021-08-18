# MAGI usage in MATLAB

Three complete, self-contained dynamic system examples are provided in the `examples/` directory.

  * Hes1 oscillation system
  * FitzHugh-Nagumo (FN) system
  * Protein transduction (PTrans) system

Each example lays out how to prepare the input to `MagiSolver`, and then use the output to plot inferred trajectories and compute parameter estimates.

Further details on how to set up your own ODE system for inference are provided following the installation instructions below, using the FitzHugh-Nagumo system  as an illustrative example.

## Installation

Using MAGI in Matlab requires the MAGI C++ library (`libcmagi.so` on Linux, `libcmagi.dll` on Windows).

We provide pre-built binaries for Windows and an automatic build script for Linux-based installs using CMake/g++.

### Linux

Compiling the `libcmagi.so` library on Linux can be done automatically using `./build.sh` in the base directory of Magi repository.

Compiling the Mex file `solveMagi.mexa64` can be done automatically using `./matab_build.sh` in this directory.

Before running `MagiSolver`, ensure `libcmagi.so` and `solveMagi.mexa64` are present in this matlabmagi directory (or another directory that is within Matlab's path).

### Windows

#### Using pre-built binaries
- For ease of use, we recommend using the pre-built binaries `libcmagi.dll` and `solveMagi.mexaw64` provided in the `windows/` subdirectory.
- The required BLAS and LAPACK libraries `libblas.dll`, `liblapack.dll` can be obtained from https://icl.cs.utk.edu/lapack-for-windows/lapack/ and downloading the `x64_dll` links.
- The GNU runtime libraries `libgcc_s_seh_64-1.dll` and `libgfortran_64-3.dll` may also be required. Copies are also provided in the `windows/` subdirectory.
- Before running `MagiSolver`, ensure `libblas.dll`, `liblapack.dll`, `libcmagi.dll`, `solveMagi.mexw64` are present in this matlabmagi folder (or another folder that is within Matlab's path). Place `libgcc_s_seh_64-1.dll` and `libgfortran_64-3.dll` in this directory as well if required. Note that if any DLLs are missing, Matlab will give the message "The specified module could not be found."

#### Manual build instructions

If you must rebuild `libcmagi.dll`, this can be done as follows:

- Install CMake from https://cmake.org/
- Install a compatible compiler, recommended to use MinGW64 x86_64-6.3.0-posix-seh-rt, obtained from https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe
- Obtain these dependent libraries: armadillo, Eigen, cppoptlib (https://github.com/PatWie/CppNumericalSolvers/tree/master), boost. All of these are header only and do not require compilation, so simply extract and place in MinGW64's `include/`.
- Get `libblas.dll`, `liblapack.dll` and place in MinGW64's `lib/`.  See https://icl.cs.utk.edu/lapack-for-windows/lapack/ for `x64_dll` download links.
- In a command prompt, go to the base directory of this Magi repository and run  `cmake -G "MinGW Makefiles"`
- Edit `cmagi\CMakeFiles\cmagi.dir\linklibs.rsp` to link to -lblas -llapack if CMake did not automatically detect them (if -larmadillo, -lpthread are not found, those flags are not needed and can be removed).
- In a command prompt, go into `cmagi/` and run `mingw32-make`.

If you then need to rebuild `solveMagi.mexaw64`, this can be done as follows:

- Ensure MinGW64 compiler support is installed (see https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)
- Start MATLAB in the matlabmagi working directory.
- Build with `mex("-IC:\MinGW\x86_64-w64-mingw32\include", "-I../include/","-L../cmagi/", "-lcmagi.dll", "src/solveMagi.cpp")`, remembering to replace the MinGW include with the appropriate directory on your system.

### Preparing the system of differential equations

First, we create a function that defines the ODEs.  It takes as input a parameter vector `theta` and matrix `x` where each column of `x` is a system component: 

```
function result = fnmodelODE(theta,x,t)
  V = x(:,1);
  R = x(:,2);
  
  result = zeros( size(x,1), size(x,2));
  result(:,1) = theta(3) * (V - V.^3 / 3.0 + R);
  result(:,2) = -1.0/theta(3) * ( V - theta(1) + theta(2) * R);
  
end
```

Next, we create functions that calculate the gradient of the ODEs, (1) with respect to the components `x`, and (2) with respect to the parameters `theta`, that provide output as 3-D arrays as follows:

```
function resultDx = fnmodelDx(theta,x,t) 
  resultDx = zeros(size(x,1),size(x,2),size(x,2));

  V = x(:,1);

  resultDx(:,1,1) = theta(3) * (1 - V.^2);
  resultDx(:,2,1) = theta(3);

  resultDx(:,1,2) = (-1.0 / theta(3));
  resultDx(:,2,2) = ( -1.0*theta(2)/theta(3) );

end

function resultDtheta = fnmodelDtheta(theta,x,t) 
  resultDtheta = zeros(size(x,1),length(theta),size(x,2));

  V = x(:,1);
  R = x(:,2);

  resultDtheta(:,3,1) = V - V.^3 / 3.0 + R;

  resultDtheta(:,1,2) =  1.0 / theta(3);

  resultDtheta(:,2,2) = -R / theta(3);
  resultDtheta(:,3,2) = 1.0/(theta(3)^2) * ( V - theta(1) + theta(2) * R);

end
```
Now we are ready to create a `struct` representing the dynamic system model to pass to `MagiSolver`.  Supply the three arguments -- ODEs, gradient with respect to `x`, gradient with respect to `theta` -- using the three functions created above.  Lastly, supply two vectors, `thetaLowerBound` and `thetaUpperBound`, that specify the lower and upper bounds on the parameters `theta`.

```
fnmodel.fOde = @fnmodelODE;
fnmodel.fOdeDx = @fnmodelDx;
fnmodel.fOdeDtheta= @fnmodelDtheta;
fnmodel.thetaLowerBound= [0 0 0];
fnmodel.thetaUpperBound= [Inf Inf Inf];
```








