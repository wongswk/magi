# do not compile wrapper.cpp from shell -- it is R wrapper
g++ hmc.cpp -o hmc.o -O2 -larmadillo
./hmc.o