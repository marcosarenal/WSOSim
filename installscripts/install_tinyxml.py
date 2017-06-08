#
# Author: Joris De Ridder
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_tinyxml.py
#
# Prerequisite: premake4
#



import os,shutil,subprocess


# Specify the dependency package name

packageName = "tinyxml"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/tinyxml/" 
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Install tinyxml

installProcedure = "cd dependencyDownloads/;   \
					mkdir tinyxml; \
					unzip tinyxml.zip  ; \
					cd %s;                      \
                    premake4 gmake --unicode --dynamic-runtime --ticpp-shared; \
                    make CONFIG=Debug;               \
                    cp -r ./lib %s;                    \
                    mkdir %s/include;                  \
                    cp *.h %s/include" %               \
                    (buildDir, installDir, installDir, installDir)

subprocess.call(installProcedure, shell=True)

      
