///////////////////////////////////////////////////////////
//  StepPhotonNoise.cpp
//  Implementation of the Class StepPhotonNoise
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



#include "StepPhotonNoise.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepPhotonNoise::StepPhotonNoise(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepPhotonNoise::~StepPhotonNoise(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the Photon Noise effect when required.
 */
void StepPhotonNoise::StepPhotonNoiseapplication(DataSet &m_DataSet)
{
            
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	
    
    //Retrieving parameters from DataSet
    numSmearingOverscanRows = p_DataSet->datasetGetnumSmearingOverscanRows(); 

            
    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();
    
    //Retrieving the smearing map from DataSet
    smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
    smearingMap = 0.0;
    smearingMap = m_DataSet.datasetGetsmearingMap();

    //Adding photon noise to each pixelMap image with a different seedRNG
    unsigned int i = 0;
    for (int x = 0; x < pixelMap.rows(); x++)
    {
        for (int y = 0; y < pixelMap.cols(); y++)
        {
            pixelMap(x, y) = Statistics::getPoisson(pixelMap(x, y), DataSet::seedRNG + (i++));
        }
    }

    //Adding photon noise to each smearing register image with a different seedRNG 
    for (int x = 0; x < pixelMap.rows(); x++)
    {
        for (int y = 0; y < numSmearingOverscanRows; y++)
        {
            smearingMap(x, y) = Statistics::getPoisson(smearingMap(x, y), DataSet::seedRNG + (i++));
        }
    }
    
            
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);


    //Setting the new calculated smearingMap into the DataSet
    m_DataSet.datasetSetsmearingMap(smearingMap);

        
        
    LogManager::log << "    Successfully added Photon Noise effect."; 
    GlobalVariables::logManager.LogManagerShowLog();     


}
//==============================================================================



//==============================================================================
/**
 * This function applies the Photon Noise effect to the WUVS processing when required.
 * The difference to the StepPhotonNoiseapplication method is that here is not applied
 * the photon noise to the smearing map.
 */
void StepPhotonNoise::StepPhotonNoiseWUVSapplication(DataSet &m_DataSet)
{
            
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	
    
               
    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();
   


    //Adding photon noise to each pixelMap image with a different seedRNG
    unsigned int i = 0;
    
    for (int x = 0; x < pixelMap.rows(); x++)
    {
        for (int y = 0; y < pixelMap.cols(); y++)
        {
            pixelMap(x, y) = Statistics::getPoisson(pixelMap(x, y), DataSet::seedRNG + (i++));
        }
    }

    
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);
      
   
    LogManager::log << "    Successfully added Photon Noise effect."; 
    GlobalVariables::logManager.LogManagerShowLog();     


}
//==============================================================================
