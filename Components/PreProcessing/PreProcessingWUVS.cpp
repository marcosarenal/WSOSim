///////////////////////////////////////////////////////////
//  PreProcessingWUVS.cpp
//  Implementation of the Class PreProcessingWUVS
//  Created on:      11-June-2015 1:59:59 PM
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



#include "PreProcessingWUVS.h"




//==============================================================================
/**
 * Constructor method
 */

PreProcessingWUVS::PreProcessingWUVS(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    
    // Initialization of the SeedRNG.
    DataSet::seedRNG = (unsigned int) time(0);
    LogManager::log << "    RNG seed value: " << DataSet::seedRNG;
    GlobalVariables::logManager.LogManagerAppendLog();
    
    // Get the predefined position from the input parameters file.
    ccdPredefinedPosition = p_DataSet->datasetGetccdPredefinedPosition();
    
    // Get the convolution method from the input parameters file.
    convolutionMethod = p_DataSet->datasetGetconvolutionMethod();
        
    // Get the output path to retrieve it
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();

    // Convert all strings to upper case
    transform(ccdPredefinedPosition.begin(), ccdPredefinedPosition.end(), ccdPredefinedPosition.begin(), ::toupper);
    transform(convolutionMethod.begin(), convolutionMethod.end(), convolutionMethod.begin(), ::toupper);
    
    // Setting the CCD positions according to the predefined options.
    PreProcessingWUVS::preProcessingCCDGetPredefinedCCDOffsets(p_DataSet);
    if (ccdPredefinedPosition == "A" || ccdPredefinedPosition == "B" || ccdPredefinedPosition == "C" || ccdPredefinedPosition == "D")
    {
        LogManager::log<< "    Using pre-defined CCD position (Only for PLATO CCD simulations)" << ccdPredefinedPosition << ".";
        GlobalVariables::logManager.LogManagerShowLog();
    }
    
    //Perform and set in the DataSet most of the CCD calculations that will be required in the processing steps.
    ParamsWUVS m_ParamsWUVS;
    m_ParamsWUVS.ParamsWUVScalculation(m_DataSet);
    
    //Calculates and sets the Background if needed.
    ParamsBackground m_ParamsBackground;
    m_ParamsBackground.ParamsBackgroundsetting(m_DataSet);
    
    //Take the corresponding FITS file
    outputPath = p_DataSet->datasetGetOutputPath();

    std::string inputFitsFile =  "/hd1/Projects/WSOSim/inputfiles/VUVES_LyA.fits";
    //std::string inputFitsFile =  "/home/pablo/Projects/Simulator/WUVS_Simulator/DATA/input_data/pablomarcosa44946/ocb6n0030_raw.fits";

    //Check that it is correctly opened
    if (FileUtilities::fileExists(inputFitsFile))
    {      
        //Read the input FITS Image from file, add the [SCI] extension and set it in the inputImage array
        FileUtilities::FileUtilitiesReadFits(inputFitsFile, inputImage); 
        //FileUtilities::FileUtilitiesReadExternalFits(inputFitsFile + "[SCI]" , inputImage); 
        FileUtilities::FileUtilitiesWriteFits(outputPath + "/" + prefix + "/" + prefix + "inputFitsFile", inputImage);
    }
        
    //Set the empty pixel map and subpixel map.
    m_DataSet.datasetSetinitPixelMap(inputImage);
    m_DataSet.datasetSetinitSubPixelMap(inputImage);        
          
      
    //Does not need to read the starfield catalogue from file.
//    ParamsStarfield m_ParamsStarfield;
//    m_ParamsStarfield.ParamsStarfieldSetting(m_DataSet);
    
    //Check if there must be performed the JITTER to set
    if (m_DataSet.datasetGetuseJitter())
    {    
        ParamsJitter m_ParamsJitter;
        m_ParamsJitter.ParamsJittersetParameters(m_DataSet);
    }
    
        
    //Calculate the CTE parameters to be set in the DataSet.
    ParamsCTE m_ParamsCTE;
    m_ParamsCTE.paramsCTEcreateCTEMap(m_DataSet);

        
    //Calculates and sets the FlatField.
    ParamsFlatField m_ParamsFlatField;
    m_ParamsFlatField.paramsFlatFieldcreateFFMap(m_DataSet);
    

}


//==============================================================================
/**
 * Destructor method
 */
PreProcessingWUVS::~PreProcessingWUVS(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 *This function provides the predefined offsets for CCD.
 */
void PreProcessingWUVS::preProcessingCCDGetPredefinedCCDOffsets(DataSet *p_DataSet)
{
	//pre-defined positions of the CCD (see manual for details)
	double ccdSizePLATO = 4510;
	double ccdSizePLATOFast = 2255;
	double pixelSizePLATO = 18;
	double pixelSizePLATOmm = pixelSizePLATO / 1000.;
	double predefinedX = 1 + pixelSizePLATOmm * ccdSizePLATO; //1mm offset of the PLATO CCDs from the optical axis, 81.18=18microns * 4510 pixels
	double predefinedXFast = 1 + pixelSizePLATOmm * ccdSizePLATOFast; //1mm offset of the PLATO CCDs from the optical axis, 18microns * 2255 pixels
	double predefinedY = 1; //1mm offset


	if (ccdPredefinedPosition == "A")
	{
		originOffsetXmm += -predefinedY;
		originOffsetYmm += predefinedX;
		ccdOrientation += 180;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATO;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "B")
	{
		originOffsetXmm += predefinedX;
		originOffsetYmm += predefinedY;
		ccdOrientation += 90;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATO;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "C")
	{
		originOffsetXmm += -predefinedX;
		originOffsetYmm += -predefinedY;
		ccdOrientation += 270;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATO;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "D")
	{
		originOffsetXmm += predefinedY;
		originOffsetYmm += -predefinedX;
		ccdOrientation += 0;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATO;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "AF")
	{
		originOffsetXmm += -predefinedY;
		originOffsetYmm += predefinedXFast;
		ccdOrientation += 180;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATOFast;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "BF")
	{
		originOffsetXmm += predefinedXFast;
		originOffsetYmm += predefinedY;
		ccdOrientation += 90;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATOFast;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "CF")
	{
		originOffsetXmm += -predefinedXFast;
		originOffsetYmm += -predefinedY;
		ccdOrientation += 270;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATOFast;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}

	if (ccdPredefinedPosition == "DF")
	{
		originOffsetXmm += predefinedY;
		originOffsetYmm += -predefinedXFast;
		ccdOrientation += 0;
		ccdSizeX = ccdSizePLATO;
		ccdSizeY = ccdSizePLATOFast;
		pixelSize = pixelSizePLATO;
                
                // If a predefined position is selected, those parameters are set into the Dataset.
                p_DataSet->datasetSetPredefinedCCDParams(originOffsetXmm,originOffsetYmm, ccdOrientation, ccdSizeX, ccdSizeY, pixelSize);
	}
}
//==============================================================================

