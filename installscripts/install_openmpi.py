#
# Author: Joris De Ridder
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_openmpi.py


import os,shutil,subprocess


# Specify the dependency package name

packageName = "openmpi-1.6.5"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/" + packageName
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Install OpenMPI

configOptions = "--disable-mpi-f77 --disable-mpi-f90"
installProcedure = "cd dependencyDownloads/;   \
					tar -xzvf openmpi-1.6.5.tar.gz; \
					cd %s;                      \
                    ./configure --prefix=%s %s; \
                    make;                       \
                    make install" % (buildDir, installDir, configOptions)

subprocess.call(installProcedure, shell=True)
                                  
