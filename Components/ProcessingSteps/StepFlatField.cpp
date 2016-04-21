///////////////////////////////////////////////////////////
//  StepFlatField.cpp
//  Implementation of the Class StepFlatField
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



#include "StepFlatField.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepFlatField::StepFlatField(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepFlatField::~StepFlatField(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function adds the FlatField effect to the image at subpixel level.
 * The flatfield map is multiplied with the sub-pixel flatfield map.
 */
void StepFlatField::StepFlatFieldapplicationSubpixel(DataSet &m_DataSet)
{
        
        //Initializing and retrieving the subpixel map from DataSet
        subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
        subPixelMap = 0.0;
        subPixelMap = m_DataSet.datasetGetsubPixelMap();
        
        //Retrieving the flatfield subpixel map from DataSet        
//        m_DataSet.datasetGetsubFlatFieldMap(subflatfieldMap);
        subflatfieldMap.resize(subPixelMap.rows(), subPixelMap.cols());
        subflatfieldMap = 0.0;
            
        subflatfieldMap = m_DataSet.datasetGetFlatFieldMap();
        

        //Application of the flatfield effect to the sub pixel map
        subPixelMap *= subflatfieldMap;
        
                
        //Setting the new calculated subPixelMap into the DataSet
        m_DataSet.datasetSetsubPixelMap(subPixelMap);
        
              
        LogManager::log <<"    Successfully added Flat Field effect at subpixel level.";
        GlobalVariables::logManager.LogManagerShowLog();     

}
//==============================================================================


//==============================================================================
/**
 * This function adds the FlatField effect to the image at pixel level.
 * The  pixelMap map is multiplied with the normal flatfield pixel map
 */
void StepFlatField::StepFlatFieldapplicationPixel(DataSet &m_DataSet)
{
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        pixelMap = m_DataSet.datasetGetpixelMap();

        //Retrieving the flatfield map from DataSet      
        flatfieldMap.resize(pixelMap.rows(), pixelMap.cols());
        flatfieldMap = 0.0;
        flatfieldMap = m_DataSet.datasetGetFlatFieldMap();
        
        
        //Application of the flatfield effect to the pixel map
        pixelMap *= flatfieldMap;

        
        //Setting the new calculated pixelMap into the DataSet
        m_DataSet.datasetSetpixelMap(pixelMap);
        
        
        LogManager::log <<"    Successfully added Flat Field effect.";
        GlobalVariables::logManager.LogManagerShowLog();     

}
//==============================================================================
