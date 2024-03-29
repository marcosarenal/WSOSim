#
# Author: Joris De Ridder
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_CCfits.py
#

import os,shutil,subprocess

# Specify the dependency package name

packageName = "CCfits"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir = currentWorkingDir + "/dependencyDownloads/" + packageName
installDir = currentWorkingDir + "/dependencyInstalls/" + packageName 
cfitsio_installDir = currentWorkingDir + "/dependencyInstalls/cfitsio"

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Install CCfits

configOptions = "--with-cfitsio=%s" % cfitsio_installDir
installProcedure = "cd dependencyDownloads/;   \
					tar -xzvf CCfits*.tar.gz; \
					cd %s;                      \
                    ./configure --prefix=%s %s; \
                    make;                       \
                    make install" % (buildDir, installDir, configOptions)

subprocess.call(installProcedure, shell=True)
            
        
