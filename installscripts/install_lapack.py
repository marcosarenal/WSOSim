#
# Author: Joris De Ridder
#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_lapack.py
#
# Prerequisite: cmake


import os,shutil,subprocess


# Specify the dependency package name

packageName = "lapack-3.5.0"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/" + packageName
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Install lapack

cmakeOptions = ""
installProcedure = "cd dependencyDownloads/;   \
					tar -xzvf lapack*.tgz; \
					cd %s;                      \
                    cmake -D CMAKE_INSTALL_PREFIX=%s %s .; \
                    make;                       \
                    make install" % (buildDir, installDir, cmakeOptions)

subprocess.call(installProcedure, shell=True)
