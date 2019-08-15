#! /bin/bash
# if the following commands error, cause the script to error
set -e

if [[ -z "${R_LIBS_USER}" ]]; then
    export R_LIBS_USER=$HOME/R/library
fi

if [[ -z "${PROJECT}" ]]; then
    PROJECT="$(git rev-parse --show-toplevel)"
    export PROJECT
fi

if [[ -z "${CPU}" ]]; then
    export CPU=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
fi

cd "$PROJECT"/rgpds || exit 1

echo "
PKG_CXX=clang++
PKG_CXXFLAGS = -std=c++11 -O3 -DNDEBUG -Wall -I$PROJECT/include
PKG_LIBS = \$(LAPACK_LIBS) \$(BLAS_LIBS) \$(FLIBS)
MAKEFLAGS = -j$CPU
" > src/Makevars

MAKE="make -j$CPU" Rscript -e 'if (!require("devtools")) install.packages("devtools", repos="http://cran.us.r-project.org")'
MAKE="make -j$CPU" Rscript -e 'if (!require("roxygen2")) install.packages("roxygen2", repos="http://cran.us.r-project.org")'
MAKE="make -j$CPU" Rscript -e 'devtools::install_deps(dependencies = TRUE, repos="http://cran.us.r-project.org")'

find src/ -lname '*' -delete
mkdir src/rcppgpds
rsync -az ../gpds_cpp/*.cpp src/rcppgpds/
rsync -az ../gpds_cpp/*.h src/rcppgpds/
perl -pi -e 's/\#include <armadillo>/\#include \"RcppArmadillo.h\"/g' src/rcppgpds/*.h
perl -pi -e 's/\#include <armadillo>/\#include \"RcppArmadillo.h\"/g' src/rcppgpds/*.cpp

perl -pi -e 's/std::cout/Rcpp::Rcout/g' src/rcppgpds/paralleltempering.cpp

rsync -az src/rcppgpds/*.cpp src/
rsync -az src/rcppgpds/*.h src/

rm -r src/rcppgpds/

Rscript -e "devtools::document();"
Rscript -e "pkgbuild:::compile_rcpp_attributes(); Rcpp::compileAttributes(); devtools::document(); devtools::install();"

LIB_PYGPDS_SOURCE=$(cd "$PROJECT"/gpds_cpp && ls -- *.cpp)
LIB_PYGPDS_HEADERS=$(cd "$PROJECT"/gpds_cpp && ls -- *.h)
cd "$PROJECT"/rgpds/src && rm $LIB_PYGPDS_SOURCE $LIB_PYGPDS_HEADERS
cd "$PROJECT"/rgpds || return
ln -s "$(pwd)"/../gpds_cpp/*.cpp src/
ln -s "$(pwd)"/../gpds_cpp/*.h src/
