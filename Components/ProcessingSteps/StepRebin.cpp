///////////////////////////////////////////////////////////
//  StepRebin.cpp
//  Implementation of the Class StepRebin
//  Created on:      23-Oct-2012 2:00:01 PM
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



#include "StepRebin.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepRebin::StepRebin(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepRebin::~StepRebin(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the rebinning to subPixelMap.
 */
void StepRebin::StepRebinapplication(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	
    
    //Retrieving parameters from DataSet
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();

    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();        

    
    //Initializing and retrieving the subpixel map from DataSet
    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
    subPixelMap = 0.0;
    subPixelMap = m_DataSet.datasetGetsubPixelMap();

            
    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();

    

    if (subPixelMap.rows() == 0 || pixelMap.rows() == 0)
    {
        return;
    }

	for (int i = 0; i < subFieldSizeX; i++)
    {
        for (int j = 0; j < subFieldSizeY; j++)
        {
            pixelMap(i, j) = sum(subPixelMap(Range(i * subPixelsPerPixel, (i + 1) * subPixelsPerPixel - 1), 
                                             Range(j * subPixelsPerPixel, (j + 1) * subPixelsPerPixel - 1)));
        }
    }
    
            
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);
    
    //Setting the new calculated subPixelMap into the DataSet
    m_DataSet.datasetSetsubPixelMap(subPixelMap);  
    
    
    LogManager::log <<"    Successfully applied rebinning.";
    GlobalVariables::logManager.LogManagerShowLog();     

}
//==============================================================================



