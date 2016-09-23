///////////////////////////////////////////////////////////
//  PreProcessingCCD.cpp
//  Implementation of the Class PreProcessingCCD
//  Created on:      23-Oct-2012 1:59:59 PM
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



#include "PreProcessingCCD.h"




//==============================================================================
/**
 * Constructor method
 */

PreProcessingCCD::PreProcessingCCD(DataSet &m_DataSet)
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
    
    // Convert all strings to upper case
    transform(ccdPredefinedPosition.begin(), ccdPredefinedPosition.end(), ccdPredefinedPosition.begin(), ::toupper);
    transform(convolutionMethod.begin(), convolutionMethod.end(), convolutionMethod.begin(), ::toupper);
    
    // Setting the CCD positions according to the predefined options.
    PreProcessingCCD::preProcessingCCDGetPredefinedCCDOffsets(p_DataSet);
    if (ccdPredefinedPosition == "A" || ccdPredefinedPosition == "B" || ccdPredefinedPosition == "C" || ccdPredefinedPosition == "D")
    {
        LogManager::log<< "    Using pre-defined CCD position " << ccdPredefinedPosition << ".";
        GlobalVariables::logManager.LogManagerShowLog();
    }
    
    //Perform and set in the DataSet most of the CCD calculations that will be required in the processing steps.
    ParamsCCD m_ParamsCCD;
    m_ParamsCCD.paramsCCDcalculation(m_DataSet);
    
    //Calculates and sets the Background if needed.
    ParamsBackground m_ParamsBackground;
    m_ParamsBackground.ParamsBackgroundsetting(m_DataSet);
    
    
    //Reads the starfield catalogue from file.
    ParamsStarfield m_ParamsStarfield;
    m_ParamsStarfield.ParamsStarfieldSetting(m_DataSet);
    
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

    //Chech if an exoplanet transit must be simulated
    performExoTransit = p_DataSet->datasetGetperformExoTransit();
    if (performExoTransit)
    {
        //Calculates and sets the Exoplanet transit parameters.
        ParamsTransit m_ParamsTransit;
        m_ParamsTransit.ParamsTransitsetting(m_DataSet);
   
    }



}


//==============================================================================
/**
 * Destructor method
 */
PreProcessingCCD::~PreProcessingCCD(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 *This function provides the predefined offsets for CCD.
 */
void PreProcessingCCD::preProcessingCCDGetPredefinedCCDOffsets(DataSet *p_DataSet)
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

