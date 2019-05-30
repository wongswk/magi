# do not compile wrapper.cpp from shell -- it is R wrapper
g++ -std=c++11 hmc.cpp -o hmc.o -O2 -larmadillo -L../../lib -I../../include/
./hmc.o

g++ -std=c++11 tgtdistr.cpp -o tgtdistr.o -O2 -larmadillo

g++ -std=c++11 hmc.cpp -o hmc.o -O2 -larmadillo -L../../lib -I../../include/

g++ -std=c++11 -c -Wall -Werror -fpic hmc.cpp -o hmc.o -O2 -I../../include/
g++ -shared -o libhmc.so hmc.o
g++ -Wall -o test hmc_main.cpp -lhmc -L../../lib -L./

g++ -std=c++11 -c -Wall -fpic tgtdistr.cpp -o tgtdistr.o -O2 -I../../include/
g++ -std=c++11 -c -Wall -fpic band.cpp -o band.o -O2 -I../../include/
g++ -shared -o libtgtdistr.so band.o tgtdistr.o -larmadillo -L../../lib -I../../include/
