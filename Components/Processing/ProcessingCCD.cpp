///////////////////////////////////////////////////////////
//  ProcessingCCD.cpp
//  Implementation of the Class ProcessingCCD
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



#include "ProcessingCCD.h"


//==============================================================================
/**
 * Constructor method
 */
ProcessingCCD::ProcessingCCD(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ProcessingCCD::~ProcessingCCD(){}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================
/**
 *This function states the processing steps to be performed in the global pipeline.
 */
void ProcessingCCD::processingCCDPipeline(DataSet &m_DataSet)
{
    
    
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving the number of exposures from the DataSet.
    numExposures = m_DataSet.datasetGetnumExposures();
    outputPath = m_DataSet.datasetGetOutputPath();
    prefix = m_DataSet.datasetGetPrefix();
    outputDir = outputPath + "/" + prefix;
    
    flatfieldIntraPixelWidth = m_DataSet.datasetGetflatfieldIntraPixelWidth();
    edgePixels = p_DataSet->datasetGetedgePixels();
    numSmearingOverscanRows  = p_DataSet->datasetGetnumSmearingOverscanRows();
    numPrescanRows = p_DataSet->datasetGetnumPrescanRows();    

    exposureTime = m_DataSet.datasetGetexposureTime();
    readOutTime = m_DataSet.datasetGetreadOutTime();
    integrationTime = exposureTime + readOutTime;    

    //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
    DataSet::seedRNG++;
    //DataSet::seedRNG  += (exposure0 * 6);
    
    //Create an array with the exposure files names
    char tempChar[18];
    exposuresNamesArray.resize(numExposures);
    for (int iter = 0; iter < numExposures; iter++)
    {
        sprintf(tempChar,  "Exposure%06d", iter);
        exposuresNamesArray(iter) = tempChar;
    }
    
    
    
    //Setting the generated FlatField map in the Dataset
    m_DataSet.datasetSetExposuresNamesArray(exposuresNamesArray);
    
    //Iterating the CCD processing for each exposure.
    for (int iterExposure = 0; iterExposure < numExposures; iterExposure++)
    {
        
        LogManager::log << "  Exposure No." << iterExposure<<": ";
        GlobalVariables::logManager.LogManagerShowLog();
        
        //Starting time for each exposure to be applied the Jitter
        startTime = iterExposure * integrationTime;
        //TODO: include exposure0 when parallel processing. exposure0 to be the first exposure in each parallel processing (pma)
        //startTime = (iterExposure + exposure0) * integrationTime;

             

        //Initialize pixelMap and subPixelMap for each new exposure
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
        subPixelMap = 0.0;
        
        //Set the null values for pixelMap and subPixelMap at the Dataset
        m_DataSet.datasetSetpixelMap(pixelMap);
        m_DataSet.datasetSetsubPixelMap(subPixelMap);
        
        
        //EXOPLANETS TRANSITS
        //Check if there TRANSITS must be used
        if (m_DataSet.datasetGetperformExoTransit())
        {    
            m_StepTransit.StepTransitapplication(m_DataSet, startTime);

        }
       
            
        //JITTER
        //Check if there JITTER must be used
        if (m_DataSet.datasetGetuseJitter())
        {         
            m_StepJitter.StepJitterapplication(m_DataSet, startTime);
        }
        else
        {
            //Retrieve the initial pixelMap and subPixelMap to be set as pixel and subPixelMap
            //This is not needed with jitter as the stars position are recalculated for each jitter position
            //Initialize the initPixelMap
            initPixelMap.resize(m_DataSet.datasetGetinitPixelMap().extent(0),m_DataSet.datasetGetinitPixelMap().extent(1));
            initPixelMap = 0.0;
            //Retrieving the initPixelMap from DataSet            
            initPixelMap = m_DataSet.datasetGetinitPixelMap();
    
            //Initialize the initSubPixelMap
            initSubPixelMap.resize(m_DataSet.datasetGetinitSubPixelMap().extent(0),m_DataSet.datasetGetinitSubPixelMap().extent(1));
            initSubPixelMap = 0.0;
            //Retrieving the initSubPixelMap from DataSet            
            initSubPixelMap = m_DataSet.datasetGetinitSubPixelMap();
 
            //Set the initial pixelMap and subPixelMap as pixel and subPixelMap
            m_DataSet.datasetSetpixelMap(initPixelMap);
            m_DataSet.datasetSetsubPixelMap(initSubPixelMap);            
        }

        //CONVOLVE PSF
        m_StepConvolvePSF.StepConvolvePSFapplication(m_DataSet);
        
        
        //CHARGE TRANSFER SMEARING
        m_StepChargeTransferSmearing.StepChargeTransferSmearingapplication(m_DataSet);
        			
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;

        //COSMICS
        //Application of COSMICS hits
        m_StepCosmics.StepCosmicsapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //BACKGROUND
        m_StepBackground.StepBackgroundapplication(m_DataSet);
        
        
        //FLATFIELD
        //Checks if the FLATFIELD map should be computed at sub-pixel level.
        //if sub-pixel white noise is > 0 or the intra-pixel width is > 0:
        if (m_DataSet.datasetGetflatfieldWhiteNoise() != 0 || m_DataSet.datasetGetflatfieldIntraPixelWidth() != 0)
        {
            m_StepFlatField.StepFlatFieldapplicationSubpixel(m_DataSet);
        }
        
          
        //REBIN
        m_StepRebin.StepRebinapplication(m_DataSet);
        
        
        //FLATFIELD
        //Checks if the FLATFIELD map should be computed at pixel level (instead of sub-pixel level).
        if (m_DataSet.datasetGetflatfieldWhiteNoise() == 0 && m_DataSet.datasetGetflatfieldIntraPixelWidth() == 0)
        {
            m_StepFlatField.StepFlatFieldapplicationPixel(m_DataSet);
        }
        

        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;

        //PHOTON NOISE
        //Check whether it is required the photon noise to be applied or not.
        if (m_DataSet.datasetGetusePhotonNoise())
        {
            m_StepPhotonNoise.StepPhotonNoiseapplication(m_DataSet);
        }
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //FULL WELL SATURATION
        m_StepSaturation.StepSaturationfullWellapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //CTE (CHARGE TRANSFER EFFICIENCY)
        m_StepCTE.StepCTEapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //READ OUT NOISE
        m_StepReadOutNoise.StepReadOutNoiseapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
        //in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //GAIN
        m_StepGain.StepGainapplication(m_DataSet);
        
        
        //ELECTRONIC OFFSET
        m_StepElectronicOffset.StepElectronicOffsetapplication(m_DataSet);
        
        
        //DIGITAL SATURATION
        m_StepSaturation.StepSaturationdigitalapplication(m_DataSet);
        
        
        
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        pixelMap = m_DataSet.datasetGetpixelMap();


        //Retrieving the smearing map from DataSet        
        //Initialize smearingMap
        smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
        smearingMap = 0.0;
        smearingMap = m_DataSet.datasetGetsmearingMap();
        
        //Retrieving the BIAS register map from DataSet        
        //Initialize smearingMap
        biasRegisterMap.resize(m_DataSet.datasetGetbiasRegisterMap().extent(0), m_DataSet.datasetGetbiasRegisterMap().extent(1));
        biasRegisterMap = 0.0;
        biasRegisterMap = m_DataSet.datasetGetbiasRegisterMap();
         
        
        //Write each exposure to FITS file cutting edges and adding BIAS register and smearing maps. 
        FileUtilities::FileUtilitiesCutAndWriteFITS(exposureTime, edgePixels, numPrescanRows, numSmearingOverscanRows, biasRegisterMap,
                                                    smearingMap, pixelMap, outputDir + "/" + prefix + exposuresNamesArray(iterExposure));

    }
      
}
//==============================================================================




//==============================================================================
/**
 *This function in a copy of the regular processingCCDPipeline. It is implemented here 
 * as a shortcut to perform the photometry on-the-go
 */
void ProcessingCCD::processingCCDPipelineOnTheGoPhotom(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry)
{
    
    
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving the number of exposures from the DataSet.
    numExposures = m_DataSet.datasetGetnumExposures();
    outputPath = m_DataSet.datasetGetOutputPath();
    prefix = m_DataSet.datasetGetPrefix();
    outputDir = outputPath + "/" + prefix;
    
    flatfieldIntraPixelWidth = m_DataSet.datasetGetflatfieldIntraPixelWidth();
    edgePixels = p_DataSet->datasetGetedgePixels();
    numSmearingOverscanRows  = p_DataSet->datasetGetnumSmearingOverscanRows();
    numPrescanRows = p_DataSet->datasetGetnumPrescanRows();    

    exposureTime = m_DataSet.datasetGetexposureTime();
    readOutTime = m_DataSet.datasetGetreadOutTime();
    integrationTime = exposureTime + readOutTime;    

	//re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
	DataSet::seedRNG++;
    //DataSet::seedRNG  += (exposure0 * 6);
    
    //Create an array with the exposure files names
    char tempChar[18];
    exposuresNamesArray.resize(numExposures);
    for (int iter = 0; iter < numExposures; iter++)
    {
        sprintf(tempChar,  "Exposure%06d", iter);
        exposuresNamesArray(iter) = tempChar;
    }
    
    
    
    //Setting the generated FlatField map in the Dataset
    m_DataSet.datasetSetExposuresNamesArray(exposuresNamesArray);
    
    //Iterating the CCD processing for each exposure.
    for (int iterExposure = 0; iterExposure < numExposures; iterExposure++)
    {
        
        LogManager::log << "  Exposure No." << iterExposure<<": ";
        GlobalVariables::logManager.LogManagerShowLog();
        
        //Starting time for each exposure to be applied the Jitter
        startTime = iterExposure * integrationTime;
        //TODO: include exposure0 when parallel processing. exposure0 to be the first exposure in each parallel processing (pma)
        //startTime = (iterExposure + exposure0) * integrationTime;

             

        //Initialize pixelMap and subPixelMap for each new exposure
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
        subPixelMap = 0.0;

        m_DataSet.datasetSetpixelMap(pixelMap);
        m_DataSet.datasetSetsubPixelMap(subPixelMap);
            
        //JITTER
        //Check if there JITTER must be used
        if (m_DataSet.datasetGetuseJitter())
        {         
            m_StepJitter.StepJitterapplication(m_DataSet, startTime);
        }
        else
        {
            //Retrieve the initial pixelMap and subPixelMap to be set as pixel and subPixelMap
            //This is not needed with jitter as the stars position are recalculated for each jitter position
            //Initialize the initPixelMap
            initPixelMap.resize(m_DataSet.datasetGetinitPixelMap().extent(0),m_DataSet.datasetGetinitPixelMap().extent(1));
            initPixelMap = 0.0;
            //Retrieving the initPixelMap from DataSet            
            initPixelMap = m_DataSet.datasetGetinitPixelMap();


            //Initialize the initSubPixelMap
            initSubPixelMap.resize(m_DataSet.datasetGetinitSubPixelMap().extent(0),m_DataSet.datasetGetinitSubPixelMap().extent(1));
            initSubPixelMap = 0.0;
            //Retrieving the initSubPixelMap from DataSet            
            initSubPixelMap = m_DataSet.datasetGetinitSubPixelMap();
            

            //Set the initial pixelMap and subPixelMap as pixel and subPixelMap
            m_DataSet.datasetSetpixelMap(initPixelMap);
            m_DataSet.datasetSetsubPixelMap(initSubPixelMap);            
            
            
            
        }

        //CONVOLVE PSF
        m_StepConvolvePSF.StepConvolvePSFapplication(m_DataSet);
        
        
        //CHARGE TRANSFER SMEARING
        m_StepChargeTransferSmearing.StepChargeTransferSmearingapplication(m_DataSet);
        			
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;

        //COSMICS
        //Application of COSMICS hits
        m_StepCosmics.StepCosmicsapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //BACKGROUND
        m_StepBackground.StepBackgroundapplication(m_DataSet);
        
        
        //FLATFIELD
        //Checks if the FLATFIELD map should be computed at sub-pixel level.
        //if sub-pixel white noise is > 0 or the intra-pixel width is > 0:
        if (m_DataSet.datasetGetflatfieldWhiteNoise() != 0 || m_DataSet.datasetGetflatfieldIntraPixelWidth() != 0)
        {
            m_StepFlatField.StepFlatFieldapplicationSubpixel(m_DataSet);
        }
        
          
        //REBIN
        m_StepRebin.StepRebinapplication(m_DataSet);
        
        
        //FLATFIELD
        //Checks if the FLATFIELD map should be computed at pixel level (instead of sub-pixel level).
        if (m_DataSet.datasetGetflatfieldWhiteNoise() == 0 && m_DataSet.datasetGetflatfieldIntraPixelWidth() == 0)
        {
            m_StepFlatField.StepFlatFieldapplicationPixel(m_DataSet);
        }
        

        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;

        //PHOTON NOISE
        //Check whether it is required the photon noise to be applied or not.
        if (m_DataSet.datasetGetusePhotonNoise())
        {
            m_StepPhotonNoise.StepPhotonNoiseapplication(m_DataSet);
        }
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //FULL WELL SATURATION
        m_StepSaturation.StepSaturationfullWellapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //CTE (CHARGE TRANSFER EFFICIENCY)
        //m_StepCTE.StepCTEapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //READ OUT NOISE
        m_StepReadOutNoise.StepReadOutNoiseapplication(m_DataSet);
        
        //re-seed the RNG with an integer incremented by 1 to ensure that each image uses a different RNG seed
		//in each noise process, each random variable is seeded with this value!
        DataSet::seedRNG++;
        
        //GAIN
        m_StepGain.StepGainapplication(m_DataSet);
        
        
        //ELECTRONIC OFFSET
        m_StepElectronicOffset.StepElectronicOffsetapplication(m_DataSet);
        
        
        //DIGITAL SATURATION
        m_StepSaturation.StepSaturationdigitalapplication(m_DataSet);
        
        //PHOTOMETRY ON-THE-GO
        m_StepPhotometryOnTheGo.StepPhotometryOnTheGoPipeline(m_DataSet, datasetPhotometry);
        	

        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        pixelMap = m_DataSet.datasetGetpixelMap();

        //Retrieving the smearing map from DataSet        
        //Initialize smearingMap
        smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
        smearingMap = 0.0;
        smearingMap = m_DataSet.datasetGetsmearingMap();
        
        //Retrieving the BIAS register map from DataSet        
        //Initialize smearingMap
        biasRegisterMap.resize(m_DataSet.datasetGetbiasRegisterMap().extent(0), m_DataSet.datasetGetbiasRegisterMap().extent(1));
        biasRegisterMap = 0.0;
        biasRegisterMap = m_DataSet.datasetGetbiasRegisterMap();
        
        //WRITE FITS FILES
        //Check whether it is required to write the FITS files to disk or not.
        bool writeFitsFilesToDisk = 1;
        if (writeFitsFilesToDisk)
        {
            //Write each exposure to FITS file cutting edges and adding BIAS register and smearing maps. 
            FileUtilities::FileUtilitiesCutAndWriteFITS(exposureTime, edgePixels, numPrescanRows, numSmearingOverscanRows, biasRegisterMap,
                                                        smearingMap, pixelMap, outputDir + "/" + prefix + exposuresNamesArray(iterExposure));
        }
        
        

    }
      
}
//==============================================================================
