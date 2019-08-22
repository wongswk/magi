#!/bin/bash

export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
export MAKE="make -j $CPU"

PROJECT=$(pwd)
ARMADILLO=armadillo-code

#make sure submodules are installed
git submodule init
git submodule sync
git submodule update

cd $PROJECT

if [ ! -d "include" ]; then
    mkdir include
fi
if [ ! -d "lib" ]; then
    mkdir lib
fi
if [ ! -d "docs" ]; then
    mkdir docs
fi

if [ ! -d "package" ]; then
    mkdir package
fi

if [ ! -f "lib/libarmadillo.dylib" ] && [ ! -f "lib/libarmadillo.so" ]; then
    cd package/armadillo-code
    git checkout 8.500.x
    ./configure
    make -j $CPU
    cp -Rf include/* ../../include/
    cp libarmadillo.* ../../lib/
    cd $PROJECT
fi

if [ ! -d "include/pybind11" ]; then
    cd package/pybind11
    git checkout v2.2.2
    cd $PROJECT
    cp -r package/pybind11/include/* include/
fi

if [ ! -d "include/boost" ]; then
    cd package/
    wget https://dl.bintray.com/boostorg/release/1.70.0/source/$BOOST.tar.gz
    tar xf $BOOST.tar.gz
    rm $BOOST.tar.gz
    cd $PROJECT
    cp -r package/$BOOST/boost include/
fi

if [ ! -d "include/Eigen" ]; then
    cd package/eigen-git-mirror
    git checkout 3.3.7
    mkdir -p install_dir
    mkdir -p build_dir && cd build_dir
    cmake .. -DCMAKE_INSTALL_PREFIX=$PROJECT/package/eigen-git-mirror/install_dir
    make install
    cp -Rf ../install_dir/include/eigen3/Eigen $PROJECT/include/
    cp -Rf ../install_dir/include/eigen3/unsupported $PROJECT/include/
    cd $PROJECT
fi

if [ ! -d "include/cppoptlib" ]; then
    cd package/CppNumericalSolvers
    git checkout dfd4686ef4cde941702024a70ac2edc73d5ee88c
    cd $PROJECT
    cp -r package/CppNumericalSolvers/include/* include/
fi
