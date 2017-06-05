#!/bin/bash
# -*- ENCODING: UTF-8 -*-

cd build
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j4
./WSOSim -s ../inputfiles/ccd_parameters.xml 
cd ..
exit
