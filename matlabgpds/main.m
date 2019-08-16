mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' src/maternCovTest.cpp;
mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' COMPFLAGS='$COMPFLAGS -std=c++11' src/basic_hmcC.cpp;
maternCovTest([1,1],[1 2; 2 1])
