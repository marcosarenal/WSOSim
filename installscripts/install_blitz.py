#
# Author: Joris De Ridder
# Usage: 
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_blitz.py


import os,shutil,subprocess

# Specify the dependency package name

packageName = "blitz-0.10"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/" + packageName
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Install blitz++

configOptions = "--disable-fortran"
installProcedure = "cd dependencyDownloads/;   \
					tar -xzvf blitz*.tar.gz; \
					cd %s;                      \
                    ./configure --prefix=%s %s; \
                    make lib;                       \
                    make install" % (buildDir, installDir, configOptions)

subprocess.call(installProcedure, shell=True)
                                  

