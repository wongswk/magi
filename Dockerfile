FROM ubuntu:18.04
LABEL maintainer="shihao.yang@isye.gatech.edu"
ENV PROJECT_DIR=/usr/src/app
WORKDIR $PROJECT_DIR

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install cmake gcc-6 g++-6 python3-tk python3-pip python3-venv python3-dev ipython3 libblas-dev liblapack-dev gfortran r-base r-base-dev git -y
#RUN add-apt-repository ppa:marutter/rrutter3.5 -y && apt-get update && apt install r-api-3.5 -y

ENV mkdir -p $PROJECT_DIR/R/library

ENV PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/lib64/pkgconfig/
ENV R_LIBS_USER=$PROJECT_DIR/R/library

RUN python3 -m pip install numpy

COPY . .
RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJECT_DIR/gpds_cpp/
RUN export PYTHONPATH=$PROJECT_DIR/pygpds

#RUN ./build.sh
