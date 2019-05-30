#!/bin/bash

./tools/dependencies.sh
cd gpds_cpp
cmake .
make
