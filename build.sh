#!/bin/bash

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

PROJECT=$(pwd)

cd $PROJECT

./tools/dependencies.sh

# build cpp
cd gpds_cpp
cmake . && make -j $CPU

cd $PROJECT

# build python
cd pygpds
pip3 install numpy nose
cmake . && make -j $CPU
python3 -c "import pygpds"
nosetests

cd $PROJECT

# build R
export CODECOV_TOKEN="7b481576-694c-4591-8370-64f61df55bdc"

echo "
PKG_CXX=clang++
PKG_CXXFLAGS = -std=c++11 -O3 -DNDEBUG -Wall \$(SHLIB_OPENMP_CXXFLAGS) -I$PROJECT/include/$BOOST
PKG_LIBS = \$(SHLIB_OPENMP_CFLAGS) \$(LAPACK_LIBS) \$(BLAS_LIBS) \$(FLIBS)

" > rgpds/src/Makevars

cd rgpds
./r_buid.sh
Rscript -e 'devtools::test()'
Rscript -e 'covr::codecov(path = "rgpds")'
