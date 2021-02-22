FROM ubuntu:18.04
LABEL maintainer="shihao.yang@isye.gatech.edu"
ENV PROJECT_DIR=/usr/src/app
WORKDIR $PROJECT_DIR

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install gcc-6 g++-6 python3-tk python3-pip python3-venv python3-dev ipython3 libblas-dev liblapack-dev gfortran git wget software-properties-common libxml2-dev build-essential libssl-dev libcurl4-openssl-dev rsync tmux -y
RUN wget https://cmake.org/files/v3.18/cmake-3.18.1-Linux-x86_64.sh && sh cmake-3.18.1-Linux-x86_64.sh  --skip-license && ln -s $(pwd)/bin/cmake /usr/local/bin/cmake
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install r-base r-base-core r-recommended r-base-dev r-cran-rgl r-cran-rjags r-cran-snow r-cran-ggplot2 r-cran-igraph r-cran-lme4 r-cran-rjava r-cran-devtools r-cran-rjava -y
#RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' && apt update && apt install cmake
#RUN add-apt-repository ppa:marutter/rrutter3.5 -y && apt-get update && apt install r-api-3.5 -y

ENV mkdir -p $PROJECT_DIR/R/library

ENV PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/lib64/pkgconfig/
ENV R_LIBS_USER=$PROJECT_DIR/R/library

RUN python3 -m pip install numpy

COPY . .
RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJECT_DIR/cmagi/
RUN export PYTHONPATH=$PROJECT_DIR/pymagi

RUN ./build.sh --skip-tests
RUN pip3 install jupyter
RUN pip3 install Pygments==2.6.1
RUN pip3 install torch pandas matplotlib
RUN pip3 install --upgrade --force jupyter-console
RUN cd $PROJECT_DIR/pymagi
# run docker using command `docker run -d -p 8888:8888 shihaoyangphd/magi:dev`
# the jupyter notebook can be used to view R results, or to run python directly
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
