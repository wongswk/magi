#!/bin/bash

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

PROJECT=$(pwd)
BOOST=boost_1_65_1


if [ ! -d "include/$BOOST" ]; then

    cd include
    wget https://downloads.sourceforge.net/project/boost/boost/1.65.1/$BOOST.tar.gz
    tar xf $BOOST.tar.gz
    rm $BOOST.tar.gz
fi

cd $PROJECT

echo "
PKG_CXX=clang++
PKG_CXXFLAGS = -std=c++11 -O3 -DNDEBUG -Wall \$(SHLIB_OPENMP_CXXFLAGS) -I$PROJECT/include/$BOOST
PKG_LIBS = \$(SHLIB_OPENMP_CFLAGS) \$(LAPACK_LIBS) \$(BLAS_LIBS) \$(FLIBS)

" > gpds/src/Makevars
