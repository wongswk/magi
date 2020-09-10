#!/usr/bin/env bash

if [[ -z "${PROJECT}" ]]; then
    PROJECT="$(git rev-parse --show-toplevel)"
    export PROJECT
fi

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

cd "$PROJECT"/matlabmagi || exit 1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJECT/cmagi
matlab -nodisplay -nosplash -nodesktop -r "mex -v '-I../include' '-I../cmagi' '-L../cmagi' '-lcmagi' GCC='g++-6' src/solveMagi.cpp;exit;"
