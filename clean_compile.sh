#!/bin/bash
# -*- ENCODING: UTF-8 -*-

cd build
rm -rf *
cmake ..
make -j4
cd ..

exit