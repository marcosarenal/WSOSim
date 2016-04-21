#
# Author: Joris De Ridder
#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_OpenCV.py
#
# Prerequisite: cmake
#

import os,shutil,subprocess



# Specify the dependency package name

packageName = "OpenCV-2.4.2"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/" + packageName+ "/build/" 
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)
os.mkdir(currentWorkingDir + "/dependencyDownloads/" + packageName)
os.mkdir(buildDir)

# Install OpenCV

cmakeOptions = "-D BUILD_opencv_highgui=OFF -D BUILD_PYTHON_SUPPORT=OFF -D BUILD_NEW_PYTHON_SUPPORT=OFF  \
				-D OPENCV_LIBPATH=$installDir/lib -D OPENCV_INCLUDEPATH=$installDir/include"
installProcedure = "cd dependencyDownloads/;   \
					tar -xjvf OpenCV*.tar.bz2; \
					cd %s;                      \
                    cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=%s %s ..; \
                    make;                       \
                    make install" % (buildDir, installDir, cmakeOptions)

subprocess.call(installProcedure, shell=True)
