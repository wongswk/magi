# do not compile wrapper.cpp from shell -- it is R wrapper
g++ -std=c++11 hmc.cpp -o hmc.o -O2 -larmadillo
./hmc.o

g++ -std=c++11 tgtdistr.cpp -o tgtdistr.o -O2 -larmadillo
