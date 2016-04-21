///////////////////////////////////////////////////////////
//  ProcessingWUVS.cpp
//  Implementation of the Class ProcessingWUVS
//  Created on:      29-May-2015 13:00:00 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
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



#include "ProcessingWUVS.h"


//==============================================================================
/**
 * Constructor method
 */
ProcessingWUVS::ProcessingWUVS(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ProcessingWUVS::~ProcessingWUVS(){}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================
/**
 *This function states the processing steps to be performed in the global pipeline
 * for the WUVS processing.
 */
void ProcessingWUVS::processingWUVSPipeline(DataSet &m_DataSet)
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
   
    //Initialize initPixelMap and pixelMap  
    initPixelMap.resize(m_DataSet.datasetGetinitPixelMap().extent(0),m_DataSet.datasetGetinitPixelMap().extent(1));
    initPixelMap = 0.0;

    pixelMap.resize(m_DataSet.datasetGetinitPixelMap().extent(0),m_DataSet.datasetGetinitPixelMap().extent(1));
    pixelMap = 0.0; 
    
    
    
    //Retrieve the initial pixelMap and subPixelMap to be set as pixel and subPixelMap
    //This is not needed with jitter as the stars position are recalculated for each jitter position
    initPixelMap = m_DataSet.datasetGetinitPixelMap();

    //Setting the generated Exposures Names Array in the Dataset
    m_DataSet.datasetSetExposuresNamesArray(exposuresNamesArray);
                
    //Iterating the WUVS processing for each exposure.
    for (int iterExposure = 0; iterExposure < numExposures; iterExposure++)
    {
        
        LogManager::log << "  Exposure No." << iterExposure<<": ";
        GlobalVariables::logManager.LogManagerShowLog();
        
        //Starting time for each exposure to be applied the Jitter
        startTime = iterExposure * integrationTime;
        //TODO: include exposure0 when parallel processing. exposure0 to be the first exposure in each parallel processing (pma)
        //startTime = (iterExposure + exposure0) * integrationTime;

      
        
        // PSF convolution is not included in WUVS Simulations since the input image already includes the PSF model.
       
                
        
        //JITTER        
        //TODO: Jitter is not included in WUVS simulations until the subpixel model is implemented (pma).
        //Check if there JITTER must be used
        if (m_DataSet.datasetGetuseJitter())
        {         
            m_StepJitter.StepJitterapplication(m_DataSet, startTime);
        }
        else
        {
            //Retrieve the initial pixelMap and subPixelMap to be set as pixel and subPixelMap
            //This is not needed with jitter as the stars position are recalculated for each jitter position
            //Initialize pixelMap and subPixelMap for each new exposure
            pixelMap = initPixelMap ;

            //Set the initial pixelMap and subPixelMap as pixel and subPixelMap
            m_DataSet.datasetSetpixelMap(pixelMap);

        
            LogManager::log <<"    Jitter effect not added to the Exposure.";
            GlobalVariables::logManager.LogManagerShowLog();    
        }
        
     
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
        
        
        //TODO: include the REBIN step once the subpixel model is implemented for WUVS (pma).
        
        
        
        
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
            //Application of the WUVS photon noise (not including photonic noise to the smearing map).
            m_StepPhotonNoise.StepPhotonNoiseWUVSapplication(m_DataSet);
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

        
        //Write each exposure to FITS file without cutting edges and adding BIAS register and smearing maps. 
        FileUtilities::FileUtilitiesWriteFits( outputDir + "/" + prefix + exposuresNamesArray(iterExposure), pixelMap );
       
    }
      
            
    FileUtilities::FileUtilitiesWriteFits(outputDir + "/" + prefix  +  "inputFitsFile_2", initPixelMap);

}
//==============================================================================


