Install PLATOSim on your local machine by following instructions in the install.txt file. This will install libraries and create a binary to run PLATOSim.

The installation of the dependencies is a local installation that doesn't require root; it shouldn't interfere with system-installed packages.

Before testing the simulator, you have to modify the file ccd_parameters.xml. The parameter <OutputPath> has to be changed to the path of an existing directory on your disk (full path required!):
The output of the simulations will be written to a directory defined by <OutputPath> and <Prefix>.
Change the parameter <FFTWisdomPath> to an existing directory. The parameters <CatalogueFileName> and <PSFFileName> must point to existing files (you can take use the included files field_RA180.0_DEC-70.0_R1.0.coo and psf_plato_6000_field_0.0_final.dat for testing).

Further information about the installation, and the User Manual can be found at the PLATO Simulator web page: 
https://fys.kuleuven.be/ster/Software/PlatoSimulator/


To start the simulation, go to the build directory and type:

./PLATOSim [-c parameter file] [[-p photometry parameters file]]
where:
-c ... Run CCD simulation
-p ... Run photometry
parameter file: An XML-file with input parameters.
photometry parameter file: An XML-file with photometry input parameters.

Examples:

Run a CCD simulation:
./PLATOSim -c /home/PLATOSim/inputfiles/ccd_parameters.xml

Run a CCD simulation and make photometry:
./PLATOSim -c /home/PLATOSim/inputfiles/ccd_parameters.xml -p /home/PLATOSim/inputfiles/photometry_parameters.xml





For questions, remarks, problems please contact Pablo Marcos Arenal <pabloMarcosArenal@ster.kuleuven.be>.
