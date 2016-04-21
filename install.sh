#!/bin/bash
#
# Just run this script:
# $ ./install.sh

#1. Install dependencies
#-----------------------

# There is one python script under installscripts/ folder for each dependency package.
# Each script Unzip/Untar packages in dependencyDownloads/ and
# installs dependencies packages in dependencyInstalls/

python ./installscripts/install_blitz.py
python ./installscripts/install_lapack.py
python ./installscripts/install_cfitsio.py
python ./installscripts/install_CCfits.py 
python ./installscripts/install_fftw.py
python ./installscripts/install_openmpi.py
python ./installscripts/install_OpenCV.py
python ./installscripts/install_tinyxml.py

# If any of these packages is already installed, just comment the 
# corresponding line (please, be sure that the commented package is correctly installed in 
# your system, otherwise the installation will raise errors!!).


#2. Build and install the simulator
#----------------------------------

mkdir build
cd build
cmake ..
make -j 4 

