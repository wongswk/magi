mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' src/maternCovTest.cpp;
maternCovTest([1,1],[1 2; 2 1])
