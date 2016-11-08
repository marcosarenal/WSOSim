///////////////////////////////////////////////////////////
//  StepTransit.cpp
//  Implementation of the Class StepTransit
//  Created on:      23-Sep-2016 11:38:59 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
 * This file is part of the WSO Simulator (WSOSim).
 *
 * WSOSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * WSOSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with WSOSim.  If not, see <http://www.gnu.org/licenses/>.
 */




#include "StepTransit.h"


//==============================================================================
/**
 * Constructor method
 */
StepTransit::StepTransit(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepTransit::~StepTransit(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the modification in the flux of the source when the transit of the exoplanet takes place.
 * The reduction in flux is calculated as the fraction of the surface of the star ocultated by the surface of the planet.
 * The transit times are retrieved from the array inTransitArray with the beguining and end of each transit as in the example:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...) 
 */
void StepTransit::StepTransitapplication(DataSet &m_DataSet, double startTime)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    hostStarTransitRA = p_DataSet->datasetGethostStarTransitRA();       
    hostStarTransitDec = p_DataSet->datasetGethostStarTransitDec();       
    hostStarRadius = p_DataSet->datasetGethostStarRadius();       
    exoplanetRadius = p_DataSet->datasetGetexoplanetRadius();       
       
       
    //Retrieving the inTransitArray from the DataSet
    inTransitArray = m_DataSet.datasetGetinTransitArray();
 
    //Initializing and retrieving the pixel map from DataSet
    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
    subPixelMap = 0.0;
    subPixelMap = m_DataSet.datasetGetsubPixelMap();


    

    //Setting the new calculated subPixelMap into the DataSet
    m_DataSet.datasetSetsubPixelMap(subPixelMap);  

    
    LogManager::log <<"    Successfully added Jitter.";
    GlobalVariables::logManager.LogManagerShowLog();
}   