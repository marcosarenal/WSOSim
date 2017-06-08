///////////////////////////////////////////////////////////
//  DataSet.h
//  Implementation of the Class DataSet
//  Created on:      23-Oct-2012 1:59:58 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
 * This file is part of the PLATO Simulator (PLATOSim).
 * Copyright 2013-2014 KU Leuven, Belgium
 *
 * PLATOSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PLATOSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with PLATOSim.  If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef DATASET_H_
#define DATASET_H_

#include <string>
#include "GlobalVariables.h"
#include "LogManager.h"
#include "ticpp.h"
#include "blitz/array.h"

using namespace std;
using namespace ticpp;
using namespace blitz;

/**
 * This class contains and stores all the parameters required to perform the simulation.
 */
class DataSet
{

public:
	
        DataSet();
        virtual ~DataSet();
        
        void datasetReadParameterFile(string parameterFile);

        
        //Parameter getters:
        inline string datasetGetOutputPath() const{ return outputPath;}                         //Retrieves the output path.
        inline string datasetGetPrefix() const{ return prefix;}                                 //Retrieves the prefix of output files
        inline string datasetGetCatalogueFileName() const{ return catalogueFileName;}           //Retrieves the stellar catalogue file name(RA, DEC, Mag).
        inline int    datasetGetnumExposures() const{ return numExposures;}                     //Retrieves the total number of exposures.

        inline string datasetGetconvolutionMethod() const{ return convolutionMethod;}           //Retrieves the convolution method (FFT or Real).
        inline int    datasetGetthreads() const{ return threads;}                               //Retrieves the number of processors to use in parallel
        inline bool   datasetGetuseFFTWisdom() const{ return useFFTWisdom;}                     //Retrieves the use-of-Wisdom value.
        inline string datasetGetfftWisdomPath() const{ return fftWisdomPath;}                   //Retrieves the directory where wisdom is stored.
        inline int    datasetGetmemoryLimit() const{ return memoryLimit;}                       //Retrieves the memory limit for FFT convolution.

        inline bool   datasetGetuseJitter() const{ return useJitter;}                           //Apply jitter to time series (0=no/1=yes).
        inline bool   datasetGetuseJitterFromFile() const{ return useJitterFromFile;}           //Apply jitter from files (0=apply jitter from parameters/1=apply jitter from file).
        inline double datasetGetjitterInterval() const{ return jitterInterval;}                 //Retrieves the time interval in between jitter positions [in seconds].
        inline double datasetGetjitterRms() const{ return jitterRms;}                           //Retrieves the jitter rms [in arcsecs].
        inline double datasetGetjitterDrift() const{ return jitterDrift;}                       //Retrieves the jitter drift [in arcsecs/min].
        inline double datasetGetjitterRepointing() const{ return jitterRepointing;}             //Retrieves the jitter repointing time [in hours]. Time for the shuttle to point back to the initial position (center of the FoV).
           
        inline bool   datasetGetusePhotonNoise() const{ return usePhotonNoise;}                   //Retrieves 0=No Photon Noise, 1=Use Photon Noise.
        inline double datasetGetreadOutNoise() const{ return readOutNoise;}                       //Retrieves the readout noise [e].
        inline int    datasetGetnumSmearingOverscanRows() const{ return numSmearingOverscanRows;} //Retrieves the smearing overscan strip number of rows [pixels]
        inline int    datasetGetnumPrescanRows() const{ return numPrescanRows;}                   //Retrieves the number of pre-scan rows [pixels].


        inline string datasetGetpsfFileName() const{ return psfFileName;}                       //Retrieves the PSF file name.
        inline double datasetGetpsfGaussFWHM() const{ return psfGaussFWHM;}                     //Retrieves the width of Gaussian PSF (only applicable if no PSF file is read in) [pixels]
        inline bool   datasetGetpsfLocationDependent()const{return psfLocationDependent;}       //Retrieves 0=no location-dependent PSF / 1=location-dependent PSF
        inline string datasetGetpsfLocationFile()const{return psfLocationFile;}                 //Retrieves the file with list of PSFs and their location.
        inline double datasetGetpsfRotationAngle() const{ return psfRotationAngle;}             //Retrieves the rotation angle of PSF in degrees (counter-clockwise).
        inline int    datasetGetpsfSubPixels() const{ return psfSubPixels;}                     //Retrieves the number of sub-pixels per pixel in PSF file (PSFNumRows/SubPixelsPerPixel must be an integer value).
        inline int    datasetGetpsfNumRows() const{ return psfNumRows;}                     //Retrieves the total number of pixels of the PSF
        inline double datasetGetpsfCenterX() const{ return psfCenterX;}                         //Retrieves the X Center of the PSF in sub-pixel coordinates.
        inline double datasetGetpsfCenterY() const{ return psfCenterY;}                         //Retrieves the Y Center of the PSF in sub-pixel coordinates.
        inline int    datasetGetpsfRotation() const{ return psfRotation;}                       //Retrieves 0 = No PSF Rotation / 1 = PSF Rotation towards Optical Axis / 2 = Arbitrary Rotation by PSFRotationAngle) 
        inline double datasetGetpsfOrientation() const{ return psfOrientation;}                 //Retrieves the orientation of pre-computed PSF in degrees (x-Axis=0°, counter-clockwise)
        inline bool   datasetGetuseGauss() const{ return useGauss;}                             //Use analytical Gaussian PSF (0=no/1=yes).

        inline string datasetGetccdPredefinedPosition() const{ return ccdPredefinedPosition;}   //Retrieves the ccd Predefined Position
        inline int    datasetGetccdSizeX() const{ return ccdSizeX;}                             //Retrieves the CCD size X [pixels] (columns).
        inline int    datasetGetccdSizeY() const{ return ccdSizeY;}                             //Retrieves the CCD size Y [pixels] (rows).
        inline double datasetGetOriginOffsetXmm() const{ return originOffsetXmm;}               //Retrieves the X Offset of CCD origin from center of focal plane [mm].
        inline double datasetGetOriginOffsetYmm() const{ return originOffsetYmm;}               //Retrieves the Y Offset of CCD origin from center of focal plane [mm].
        inline double datasetGetccdOrientation() const{ return ccdOrientation;}                 //Retrieves the orientation of CCD with respect to focal plane orientation [deg].
        inline double datasetGetorientationFocalPlane() const{ return orientationFocalPlane;}   //Retrieves the orientation of focal plane [deg].
        
        inline int    datasetGetsubFieldZeroX() const{ return subFieldZeroX;}                   //Retrieves the sub-field zero point X relative to CCD [pixels].
        inline int    datasetGetsubFieldZeroY() const{ return subFieldZeroY;}                   //Retrieves the sub-field zero point Y relative to CCD [pixels].
        inline int    datasetGetsubPixelsPerPixel() const{ return subPixelsPerPixel;}           //Retrieves the number of sub-pixels per pixel.
        inline int    datasetGetsubFieldSizeX() const{ return subFieldSizeX;}                   //Retrieves the sub-field size X [pixels] (columns)-
        inline int    datasetGetsubFieldSizeY() const{ return subFieldSizeY;}                   //Retrieves the sub-field size Y [pixels] (columns)-

        inline double datasetGetraCenterSubField() const{ return raCenterSubField;}             //Retrieves the Right Ascension of center of FOV of CCD [deg].
        inline double datasetGetdeclCenterSubField() const{ return declCenterSubField;}         //Retrieves the declination of center of FOV of CCD [deg].
        inline int    datasetGetedgePixels() const{ return edgePixels;}                         //Retrieves the edgePixels parameter.

        inline int    datasetGetpixelSize() const{ return pixelSize;}                           //Retrieves the pixel size [micron].
        inline double datasetGetpixelScale() const{ return pixelScale;}                         //Retrieves the pixel scale [arcsec/pixel]
        
        inline double datasetGetopticalAxisRACenter()const{ return opticalAxisRACenter;}        //Retrieves the Right Ascension of optical axis [deg].
        inline double datasetGetopticalAxisDecCenter()const{ return opticalAxisDecCenter;}      //Retrieves the declination of optical axis [deg].
        inline double datasetGetraOpticalAxis()const{ return raOpticalAxis;}                    //Retrieves the RA orientation of optical axis after displacement (jitter) [rad]
        inline double datasetGetdecOpticalAxis()const{ return decOpticalAxis;}                  //Retrieves the declination orientation of optical axis after displacement (jitter) [rad]
        inline double datasetGetrotationAngleOA()const{ return rotationAngleOA;}                //Retrieves the rotation angle of optical axis after displacement (jitter) [rad]

        inline double datasetGetxOpticalAxis()const{ return xOpticalAxis;}                      //Retrieves the x position of the center of the focal plane.
        inline double datasetGetyOpticalAxis()const{ return yOpticalAxis;}                      //Retrieves the y position of the center of the focal plane.
 
        inline double datasetGetxFOVSubField()const{ return xFOVSubField;}                      //Retrieves the Field of view of sub-image in x-direction.
        inline double datasetGetyFOVSubField()const{ return yFOVSubField;}                      //Retrieves the Field of view of sub-image in y-direction.
        inline int    datasetGetradiusFOVCCD() const{ return radiusFOVCCD;}                     //Retrieves the radius of the FOV [deg].
        
        inline double datasetGetareaTelescope() const{ return areaTelescope;}                   //Retrieves the light collecting area of each telescope [cm^2].
        inline double datasetGettransEff() const{ return transEff;}                             //Retrieves the transmission efficiency.
        inline double datasetGetquantEff() const{ return quantEff;}                             //Retrieves the quantum efficiency.
        inline double datasetGetbackground() const{ return background;}                         //Retrieves the Sky background (zodiacal+galactic) [e/(s*pixel)].
        inline double datasetGetfluxm0() const{ return fluxm0;}                                 //Retrieves the flux of m=0 star [phot/s/cm^2]
        inline int    datasetGetfullWellSat() const{ return fullWellSat;}                       //Retrieves the full well saturation [e/pixel]
        inline int    datasetGetdigitalSat() const{ return digitalSat;}                         //Retrieves the digital saturation [ADU/pixel]
        inline double datasetGetgain() const{ return gain;}                                     //Retrieves the gain [e/ADU].
        inline double datasetGetelectronicOffset() const{ return electronicOffset;}             //Retrieves the electronic offset (= bias level) [ADU]
        
        inline double datasetGetexposureTime() const{ return exposureTime;}                     //Retrieves the exposure time [s].
        inline bool   datasetGetstackExposures() const{ return stackExposures;}                     //Retrieves the exposure time [s].
        inline int    datasetGetnumExposuresPerStack() const{ return numExposuresPerStack;}                     //Retrieves the exposure time [s].
        inline double datasetGetreadOutTime() const{ return readOutTime;}                       //Retrieves the read-out time [s].
        inline double datasetGetcosmicHitRate() const{ return cosmicHitRate;}                   //Retrieves the cosmic hit rate [events/cm^2/min].
        inline double datasetGetcosmicsWidth() const{ return cosmicsWidth;}                     //Retrieves the cosmic hit width [pixels].
        inline double datasetGetcosmicsLength() const{ return cosmicsLength;}                   //Retrieves the cosmic hit width [pixels].
        inline double datasetGetcosmicsSatFactor() const{ return cosmicsSatFactor;}             //Retrieves the cosmic saturation factor [proportional to full-well saturation].
        
        inline double datasetGetmeanCTE() const{ return meanCTE;}                               //Retrieves the mean Charge Transfer Efficiency.
        inline int    datasetGetnumLowCTEPixels() const{ return numLowCTEPixels;}               //Retrieves the number of low-CTE pixels.
        inline int    datasetGetnumLowCTELines() const{ return numLowCTELines;}                 //Retrieves the number of low-CTE lines.

        
        inline double datasetGetflatfieldPixelNoise() const{ return flatfieldPixelNoise;}       //Retrieves the Flatfield peak-to-peak pixel noise. 
        inline double datasetGetflatfieldWhiteNoise() const{ return flatfieldWhiteNoise;}       //Retrieves the Flatfield sub-pixel white noise.
        inline double datasetGetflatfieldIntraPixelWidth() const{ return flatfieldIntraPixelWidth;}//Retrieves the Flatfieldwidth of the central part of the pixel which is affected by
                                                                                                   //a loss of sensitivity  lower than 5\% due to edge effect [\% of pixel size, rounded up]\\.

        inline string datasetGetjitterFile() const{ return jitterFile;}                         //Retrieves the file name with jitter time series
        inline double datasetGetjitterAngularDist() const{ return jitterAngularDist;}           //Retrieves the angular distance [deg] of satellite jitter rotation axis from CCD center-of-field
        inline double datasetGetjitterPosAngle() const{ return jitterPosAngle;}                 //Retrieves the position angle [deg] of satellite jitter rotation axis (relative to N).
        inline double datasetGetjitterMultFactor() const{ return jitterMultFactor;}             //Retrieves the multiplication factor of jitter time-series.
        
        inline bool   datasetGetperformExoTransit() const{ return performExoTransit;}           //Retrieves 0 = Do not perform transit simulation of exoplanet, 1 = Do perform transit simulation of exoplanet.
        inline double datasetGethostStarTransitRA() const{ return hostStarTransitRA;}           //Retrieves the Right Ascension of the transit host star [deg].
        inline double datasetGethostStarTransitDec() const{ return hostStarTransitDec;}         //Retrieves the Declination of the transit host star [deg].
        inline double datasetGethostStarRadius() const{ return hostStarRadius;}                 //Radius of the transit host star [Solar radius].
        inline double datasetGetexoplanetRadius() const{ return exoplanetRadius;}               //Retrieves the Radius of the exoplanet [Solar radius].
        inline double datasetGetexoplanetOrbitalPeriod() const{ return exoplanetOrbitalPeriod;} //Retrieves the Orbital period of the exoplanet [days].
        inline double datasetGetplanetaryOrbitSemiaxis() const{ return planetaryOrbitSemiaxis;} //Retrieves the semiaxis of the orbit of the exoplanet [AU].
        inline double datasetGetplanetaryOrbitInclination() const{ return planetaryOrbitInclination;} //Retrieves the inclination of the orbit of the exoplanet as seen from Earth [deg].
        inline double datasetGethostStarMagnitude() const{ return hostStarMagnitude;}           //Retrieves the magnitude of the transit host star on the Star Catalogue
        inline double datasetGethostStarID() const{ return hostStarID;}                         //Retrieves the Position of the transit host star on the Star Catalogue

        
        Array<float, 2>  datasetGetPSFMap();                                
        Array<double, 2> datasetGetCTEMap();                               
        Array<float, 2>  datasetGetpixelMap();                            
        Array<float, 2> datasetGetsubPixelMap(); 
        Array<float, 2> datasetGetinitPixelMap();    
        Array<float, 2> datasetGetinitSubPixelMap();                      
        Array<float, 2> datasetGetFlatFieldMap();
        Array<float, 2> datasetGetsubFlatFieldMap();              
        Array<float, 2> datasetGetsmearingMap();                      
        Array<float, 2> datasetGetbiasRegisterMap();              
        Array<double, 2>  datasetGetJitterInputParams();         
        Array<float, 2>  datasetGetStarListOnCCD();                  
        Array<float, 2>  datasetGetsubPixelStarListOnCCD();                  
        Array<string, 1>  datasetGetExposuresNamesArray();        
        Array<float, 2> datasetGetStarCatalogue();                  
        Array<float, 1> datasetGetinTransitArray();                  

        
        
        

        
        //Parameter setters:
        void datasetSetPredefinedCCDParams(int originOffsetXmm, int originOffsetYmm, int ccdOrientation, 
                int ccdSizeX, int ccdSizeY, int pixelSize);
        
        void datasetSetCCDParams(int ccdSizeX, int ccdSizeY, int subFieldSizeX, int subFieldSizeY, 
                double originOffsetXmm, double originOffsetYmm, double raCenterSubField, double declCenterSubField, 
                double radiusFOVCCD, int edgePixels, double raOpticalAxis, double decOpticalAxis, double rotationAngleOA,
                double xFOVSubField, double yFOVSubField);
        
        void datasetSetbackground(float newBackground);
        void datasetSetpixelMap(Array<float, 2> pixelMap);
        void datasetSetsubPixelMap(Array<float, 2> subPixelMap);
        void datasetSetinitPixelMap(Array<float, 2> initPixelMap);
        void datasetSetinitSubPixelMap(Array<float, 2> initSubPixelMap);
        void datasetSetPSFMap(Array<float, 2> psfMap);
        void datasetSetCTEMap(Array<double, 2> cteMap);
        void datasetSetFlatFieldMap(Array<float, 2> flatfieldMap);
        void datasetSetsubFlatFieldMap(Array<float, 2> subflatfieldMap);
        void datasetSetsmearingMap(Array<float, 2> smearingMap);
        void datasetSetbiasRegisterMap(Array<float, 2> biasRegisterMap);
        void datasetSetJitterInputParams(Array<double, 2> jitterInputParams);
        void datasetSetStarListOnCCD(Array<float, 2> starListOnCCD);
        void datasetSetsubPixelStarListOnCCD(Array<float, 2> subPixelStarListOnCCD);
        void datasetSetCenterFocalPlane(double xOpticalAxis, double yOpticalAxis);
        void datasetSetExposuresNamesArray(Array<string, 1> exposuresNamesArrayCopy);
        void datasetSetStarCatalogue(Array<float, 2> starCatalogue);
        void datasetSetinTransitArray(Array<float, 1> inTransitArray);
        void datasetSethostStarMagnitude(double hostStarMagnitude);
        void datasetSethostStarID(double hostStarID);
              


        // Static parameters:

        static unsigned int seedRNG;



private:
                
        string outputPath;                      //Output path.
        string catalogueFileName;               //File name stellar catalogue (RA, DEC, Mag).
        string prefix;                          //Prefix of output files (=name of output sub-directory).
        int    threads;                         //Number of processors to use in parallel.
        int    memoryLimit;                     //Memory limit for FFT convolution (to avoid memory overflow).
        int    edgePixels;                      //Number of pixels by which the image size is increased on each side to ensure that the treatment of the edges is made correctly.
        
        string ccdPredefinedPosition;           //Predefined position of CCD (A, B, C, D, AF, BF, CF, DF, or User).
        int    ccdSizeX;                        //CCD size X [pixels] (columns)
        int    ccdSizeY;                        //CCD size Y [pixels] (rows)
        double originOffsetXmm;                 //X Offset of CCD origin from center of focal plane [mm]
        double originOffsetYmm;                 //Y Offset of CCD origin from center of focal plane [mm]
        double ccdOrientation;                  //Orientation of CCD with respect to focal plane orientation [deg].
        double opticalAxisRACenter;             //Right Ascension of optical axis [deg].
        double opticalAxisDecCenter;            //Declination of optical axis [deg].
        double orientationFocalPlane;           //Orientation of focal plane [deg].
        double raOpticalAxis;                   //Right ascension orientation of optical axis after displacement (jitter)
        double decOpticalAxis;                  //Declination orientation of optical axis after displacement (jitter)
        double rotationAngleOA;                 //Rotation angle of optical axis after displacement (jitter)
        double xOpticalAxis;                    //x position of the center of the focal plane.
        double yOpticalAxis;                    //y position of the center of the focal plane.
        
        int    subFieldZeroY;                   //Sub-field zero point X relative to CCD [pixels].
        int    subFieldZeroX;                   //Sub-field zero point Y relative to CCD [pixels].
        int    subFieldSizeX;                   //Sub-field size X [pixels] (columns).
        int    subFieldSizeY;                   //Sub-field size X [pixels] (columns).
        int    subPixelsPerPixel;               //Number of sub-pixels per pixel.
        double xFOVSubField;                    //Field of view of sub-image in x-direction.
        double yFOVSubField;                    //Field of view of sub-image in y-direction.
        double raCenterSubField;                //Right Ascension of center of FOV of CCD [deg].
        double declCenterSubField;              //Declination of center of FOV of CCD [deg].
        double radiusFOVCCD;                    //Radius of the FOV [deg].
        
        string psfFileName;                     //PSF file name.
        string psfLocationFile;                 //File with list of PSFs and their location.
        double psfGaussFWHM;                    //Width of Gaussian PSF (only applicable if no PSF file is read in) [pixels]
        double psfCenterX;                      //X Center of the PSF in sub-pixel coordinates.
        double psfCenterY;                      //X Center of the PSF in sub-pixel coordinates.
        double psfOrientation;                  //Orientation of pre-computed PSF in degrees (x-Axis=0°, counter-clockwise).
        double psfRotationAngle;                //Rotation angle of PSF in degrees (counter-clockwise).
        int    psfRotation;                     //Rotate PSF  (0 = No Rotation / 1 = Towards Optical Axis / 2 = Arbitrary Rotation by PSFRotationAngle).
        int    psfNumRows;                    //Total number of pixels of the PSF taken from the number of rows of PSF file (must be quadratic).
        int    psfSubPixels;                    //Number of sub-pixels per pixel in PSF file (PSFNumRows/SubPixelsPerPixel must be an integer value).
        bool   psfLocationDependent;            //Use location-dependent PSF (0=no/1=yes).
        
        string convolutionMethod;               //Convolution method: FFT Convolution (FFT) or Real Space Convolution (Real).
        bool   useFFTWisdom;                    //Use FFT Wisdom to speed up convolution (can take longer for certain image dimensions) (0=no/1=yes)
        string fftWisdomPath;                   //The wisdom information is stored in this directory.
        bool   useGauss;                        //Use analytical Gaussian PSF (0=no/1=yes).
        
        double exposureTime;                    //Exposure time [s].
        int    numExposures;                    //Total number of exposures.
        bool   stackExposures;                  //Stack series of exposures in one output FITS image file (0=no/1=yes).        
        int    numExposuresPerStack;            //Number of exposures to be summed in the same image.
        double fluxm0;                          //Flux of m=0 star [phot/s/cm^2].
        
        double gain;                            //Gain [e/ADU].
        double pixelScale;                      //Pixel scale [arcsec/pixel].
        double transEff;                        //Transmission efficiency.
        double quantEff;                        //Quantum efficiency.
        double areaTelescope;                   //Light collecting area of each telescope [cm^2].
        double pixelSize;                       //Pixel size [micron].
        double readOutTime;                     //Time required to read-out a entire CCD working in frame transfer mode with open shutter [s].
        double electronicOffset;                //Electronic offset (= bias level) [ADU].
        double background;                      //Sky background (zodiacal+galactic) [e/(s*pixel)].
        double flatfieldPixelNoise;             //Flatfield peak-to-peak pixel noise.
        double flatfieldWhiteNoise;             //Flatfield sub-pixel white noise.
        double flatfieldIntraPixelWidth;        //Retrieves the Flatfieldwidth of the central part of the pixel which is affected by
                                                // a loss of sensitivity  lower than 5\% due to edge effect [\% of pixel size, rounded up]\\.
        
        double cosmicHitRate;                   //Cosmic hit rate [events/cm^2/min].
        double cosmicsSatFactor;                //Cosmics saturation factor [proportional to full-well saturation].
        double cosmicsWidth;                    //Cosmics hit width FWHM [pixels].
        double cosmicsLength;                   //Cosmics hit length FWHM [pixels].
        
        int    fullWellSat;                     //Full well saturation [e/pixel].
        int    digitalSat;                      //Digital saturation [ADU/pixel].
        double readOutNoise;                    //Readout noise [e].
        int    numSmearingOverscanRows;         //Smearing overscan strip number of rows [pixels].
        int    numPrescanRows;                  //Number of pre-scan rows [pixels].
        bool   usePhotonNoise;                  //0=No Photon Noise, 1=Use Photon Noise.
        
        bool   useJitter;                       //Apply jitter to time series (0=no/1=yes).
        bool   useJitterFromFile;               //Apply jitter to time series from a file with jitter time-series or using jitter parameters (0=using parameters/1=from file).
        double jitterInterval;                  //Retrieves the time interval in between jitter positions [in seconds].
        double jitterRms;                       //Retrieves the jitter rms [in arcsecs].
        double jitterDrift;                     //Retrieves the jitter drift [in arcsecs/min].
        double jitterRepointing;                //Retrieves the jitter repointing time [in hours]. Time for the shuttle to point back to the initial position (center of the FoV).

        string jitterFile;                      //File with jitter time series.
        double jitterPosAngle;                  //Position angle [deg] of satellite jitter rotation axis (relative to N).
        double jitterAngularDist;               //Angular distance [deg] of satellite jitter rotation axis from CCD center-of-field.
        double jitterMultFactor;                //Multiplication factor of jitter time-series.
        
        double meanCTE;                         //Mean Charge Transfer Efficiency.
        int    numLowCTELines;                  //Number of low-CTE lines.
        int    numLowCTEPixels;                 //Number of low-CTE pixels.
        
        bool   writeSubPixelFits;               //Write subpixel map to FITS file (0=no/1=yes).

        bool   performExoTransit;               //Perform exoplanetary transit simulation for one source in the field (0=no/1=yes).
        double hostStarTransitRA;               //Right Ascension of the transit host star [deg] (It must match with a source in the star catalogue).
        double hostStarTransitDec;              //Declination of the transit host star [deg].
        double hostStarRadius;                  //Radius of the transit host star [Solar radius].
        double exoplanetRadius;                 //Radius of the exoplanet [Solar radius].
        double exoplanetOrbitalPeriod;          //Orbital period of the exoplanet [days].
        double planetaryOrbitSemiaxis;          //Semiaxis of the orbit of the exoplanet [AU].
        double planetaryOrbitInclination;       //Inclination of the orbit of the exoplanet as seen from Earth [deg].
        double hostStarMagnitude;               //Magnitude of the transit host star
        double hostStarID;                      //Position of the transit host star on the Star Catalogue
        
        
        
        Array<float, 2>  pixelMap;                      //Blitz 2-D array with the pixel map.
        Array<float, 2>  subPixelMap;                   //Blitz 2-D array with the subpixel map.
        Array<float, 2>  initPixelMap;                  //Blitz 2-D array with the pixel map.
        Array<float, 2>  initSubPixelMap;               //Blitz 2-D array with the subpixel map.
        Array<float, 2>  smearingMap;                   //Blitz 2-D array with the CCD shielded register containing the star smearing due to frame transfer.
        Array<float, 2>  psfMap;                        //Blitz 2-D array with the PSF map.
        Array<double, 2> cteMap;                        //Blitz 2-D array with the CTE map.
        Array<float, 2>  flatfieldMap;                  //Blitz 2-D array with the FlatField map.
        Array<float, 2>  subFlatFieldMap;               //Blitz 2-D array with the subFlatField map.
        Array<float, 2>  biasRegisterMap;               //Blitz 2-D array with the bias overscan strip of the CCD map. 
        Array<double, 2> jitterInputParams;             //Blitz 2-D array with the jitter displacement time [s], yaw, pitch and roll angles [deg].
        Array<float, 2>  starListOnCCD;                 //Blitz 2-D array with the list of stars on the CCD.
        Array<float, 2>  subPixelStarListOnCCD;         //Blitz 2-D array with the list of stars on the CCD at subpixel level.
        Array<string, 1> exposuresNamesArray;           //Blitz 2-D array with the list of exposures names.
        Array<float, 2>  starCatalogue;                 //Blitz 2-D array with the star catalogue
        Array<float, 1>  inTransitArray;                //Blitz 2-D array with the the times where the exoplanet transit takes place.
};



/**
 * This class contains and stores all the parameters required to perform the simulation.
 */
class DataSetPhotometry
{

public:
	
        DataSetPhotometry();
        virtual ~DataSetPhotometry();
        
        void datasetPhotometryReadParameterFile(string photometryParameterFile);
        
        //Photometry Parameter getters:
        inline string datasetGetphotometryMethod() const{ return photometryMethod;}   //Retrieves the photometry method selected.
        inline string datasetGetphotometryDirName() const{ return photometryDirName;}   //Retrieves the photometry directory name.
        inline string datasetGetphotometryPlotsDir() const{ return photometryPlotsDir;}   //Retrieves the photometry directory name for plots.
        inline bool datasetGetframeTransferSmearingCorrection() const{ return frameTransferSmearingCorrection;}   //Retrieves whether the frame Transfer Smearing must be corrected or not.
        inline int datasetGetphotometryNumTelescopes() const{ return photometryNumTelescopes;}   //Retrieves the photometry number of telescopes.
        inline bool datasetGetphotometryFlatfieldCorrectionr() const{ return flatfieldCorrection;}   //Retrieves whether the flatfield must be corrected or not.
        inline double datasetGetphotometryBackground() const{ return photometryBackground;}   //Retrieves the photometry Background Flux.
        inline double datasetGetphotometryBackgroundAnnulusInnerRadius() const{ return backgroundAnnulusInnerRadius;}   //Retrieves the photometry Background Annulus Inner Radius.
        inline double datasetGetphotometryBackgroundAnnulusOuterRadius() const{ return backgroundAnnulusOuterRadius;}   //Retrieves the photometry Background Annulus Outer Radius.

        //Photometry Parameter setters:
        void datasetSetPhotometryParams(string PhotomOutputDir, string PhotomPlotsDir);
    

        
private:
                
        string photometryMethod;                //Photometry method selected (WM=Weighted Mask, AP= Simple Aperture Photometry).
        int    photometryNumTelescopes;         //Number of telescopes that should be considered for statistical calculations.
        bool   flatfieldCorrection;             //Correct the flat field (0=no/1=yes) (flatfield must be provided and have mean flux of 30000 ADU)
        bool   frameTransferSmearingCorrection; //Correct for the smearing trails from the read-out with open shutter (0=no/1=yes) (trails must be provided in a 1D file).
        double photometryBackground;            //Background flux [e-/pixel/s] (negative for automatic calculation).
        double backgroundAnnulusInnerRadius;    //Inner radius of annulus for background calculation [pixels].
        double backgroundAnnulusOuterRadius;    //Outer radius of annulus for background calculation [pixels]
        string photometryDirName;               //Directory name for the photometry output.
        string photometryPlotsDir;              //Directory name for the photometry plots.


};
#endif /* DATASET_H_ */       