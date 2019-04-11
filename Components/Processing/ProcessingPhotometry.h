///////////////////////////////////////////////////////////
//  ProcessingPhotometry.h
//  Implementation of the Class ProcessingPhotometry
//  Created on:      23-Oct-2012 2:00:00 PM
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



#ifndef PROCESSINGPHOTOMETRY_H_
#define PROCESSINGPHOTOMETRY_H_


#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"
#include "MathTools.h"
#include "Statistics.h"




/**
 * Class which calls the main methods for performing the photometry processing.
 */
class ProcessingPhotometry
{

public:
    ProcessingPhotometry();
    virtual ~ProcessingPhotometry();
    void processingPhotometryPipeline(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry);
    void processingPhotometrySubPixelStarListOnCCDCut();
    void processingPhotometryPreComputeWeightedMask();
    void processingPhotometryMeasureStarsFluxes(Array<float, 2> flux, int frameIndex);
    void processingPhotometryComputeStatistics(std::string fileName);
    void processingPhotometryMakeAnalysis(std::string fileName);
    void processingPhotometryWritePsfToASCIIFile(std::string fileName);
    void processingPhotometryWriteInfo(std::string fileName);
    
    
private:
        
    DataSetPhotometry*  p_datasetPhotometry;          //Pointer to the DataSet to retrieve parameters from it.
    DataSet*    p_DataSet;                            //Pointer to the DataSet to retrieve parameters from it.

    Array<float, 2> flux;
    Array<float, 2> rawImage;
    Array<std::string, 1> exposuresNamesArray;             //Array containing the exposures names for the FITS files.
    Array<float, 2>  subPixelStarListOnCCD;
    Array<float, 2>  psfMap;                          //Blitz 2-D array with the PSF map.
    Array<float, 4>  weightedPSF;                     //Blitz 4-D array. Each pixel of the 2-D PSF map is associated to a 2-D array 
    Array<float, 2>  normalizationFactor;
    Array<float, 2>  flatfieldMap;                    //Blitz map retrieved from the DataSet.

    std::string  outputPath;                               //Parameter retrieved from DataSet.
    std::string  prefix;                                   //Parameter retrieved from DataSet.
    std::string  photometryDirName;                        //Parameter retrieved from DataSetPhotometry.
    std::string  photometryPlotsDir;                       //Parameter retrieved from DataSetPhotometry.
    std::string  photometryMethod;                         //Parameter retrieved from DataSetPhotometry.
    
    int     numSmearingOverscanRows;                  //Parameter retrieved from DataSet.
    int     numPrescanRows;                           //Parameter retrieved from DataSet.
    double  gain;                                     //Parameter retrieved from DataSet.
    int     photometryNumTelescopes;                  //Parameter retrieved from DataSetPhotometry.
    bool    frameTransferSmearingCorrection;          //Parameter retrieved from DataSetPhotometry.
    bool    flatfieldCorrection;                      //Parameter retrieved from DataSetPhotometry.
    bool    useGauss;                                 //Parameter retrieved from DataSet.
    int     numExposures;                             //Parameter retrieved from DataSet.

    int     subPixelsPerPixel;                        //Parameter retrieved from DataSet.
    int     psfSubPixels;                             //Parameter retrieved from DataSet.
    double  backgroundAnnulusInnerRadius;             //Parameter retrieved from DataSetPhotometry.
    double  backgroundAnnulusOuterRadius;             //Parameter retrieved from DataSetPhotometry.
    double  photometryBackground;                     //Parameter retrieved from DataSetPhotometry.
    int     edgePixels;                               //Parameter retrieved from DataSet.
    double  quantEff;                                 //Parameter retrieved from DataSet.
    double  transEff;                                 //Parameter retrieved from DataSet.
    double  fluxm0;                                   //Parameter retrieved from DataSet.
    double  areaTelescope;                            //Parameter retrieved from DataSet.
    double  readOutTime;                              //Parameter retrieved from DataSet.
    double  exposureTime;                             //Parameter retrieved from DataSet.
    int     subFieldSizeX;                            //Parameter retrieved from DataSet.
    int     subFieldSizeY;                            //Parameter retrieved from DataSet.
    int     digitalSat;                               //Parameter retrieved from DataSet.       
    int     psfNumRows;                             //Parameter retrieved from DataSet.
    
    double  magObs;                                   //Magnitude observed for each star.
    double  magMC;                                    //Mean input magnitude generated by Monte-Carlo Method for each star.
    double  magInput;                                 //Input magnitude for each star.
    double  magObsNorm;                               //Mean  observed normalized magnitude for each star.
    double  fluxObs;                                  //Flux observed for each star.
    double  fluxInput;                                //Input flux for each star.
    double  fluxObsNorm;                              //Mean  observed normalized flux for each star.
    
};
#endif /* PROCESSINGPHOTOMETRY_H_ */
