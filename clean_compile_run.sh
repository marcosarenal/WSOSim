#!/bin/bash
# -*- ENCODING: UTF-8 -*-

cd build
rm -rf *
cmake ..
make -j4
./PLATOSim -w ../inputfiles/WUVS_parameters.xml 
cd ..
exit