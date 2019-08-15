#!/bin/bash
# if the following commands error, cause the script to error
set -e

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

if [[ -z "${R_LIBS_USER}" ]]; then
    export R_LIBS_USER=$HOME/R/library
fi

export PROJECT=$(pwd)
export BOOST=boost_1_70_0

cd $PROJECT

./tools/dependencies.sh

# build cpp
cd gpds_cpp
cmake . && make -j $CPU

cd $PROJECT

# build python
cd pygpds
pip3 install -r pip/requirements.txt
cmake . && make -j $CPU
python3 -c "import pygpds"
nosetests

cd $PROJECT

# build R
export CODECOV_TOKEN="7b481576-694c-4591-8370-64f61df55bdc"

cd rgpds
./r_build.sh
# Rscript -e 'devtools::test()'
Rscript -e 'testthat::test_package("gpds")'
Rscript -e 'covr::codecov(path = ".")'

cd $PROJECT
