mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' src/maternCovTest.cpp;
maternCovTest([1,1],[1, 2; 2, 1])

mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/basic_hmcC.cpp;
tgt = @loglik;
[final, lpx] = basic_hmcC(tgt, [0.2, 0.2], [0.1, 0.1], [-10, -10], [inf, inf], 200, false)
assert(abs(tgt(final) - lpx) < 1e-4)
