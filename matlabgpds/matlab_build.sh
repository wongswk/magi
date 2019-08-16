#!/usr/bin/env bash

if [[ -z "${PROJECT}" ]]; then
    PROJECT="$(git rev-parse --show-toplevel)"
    export PROJECT
fi

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

cd "$PROJECT"/matlabgpds || exit 1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJECT/gpds_cpp
matlab -nodisplay -nosplash -nodesktop -r "mex -v '-I../include' '-I../gpds_cpp' '-L../gpds_cpp' '-lcgpds' GCC='g++-6' src/maternCovTest.cpp;exit;"
