#!/bin/bash
# -*- ENCODING: UTF-8 -*-

cd build
make -j4


START_TIME=$SECONDS
./PLATOSim -w ../inputfiles/VUVES_Sp-A0I.XML 
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo ""
echo "-------------->VUVES Simulation Time= "$ELAPSED_TIME " seg."
echo ""

cd ..

#run python script to check difference between input image and Exposure 0
python ../Simulator/WUVS_Simulator/WUVSevaluation/WUVS_evaluation/WUVSevaluation.py _test

exit