///////////////////////////////////////////////////////////
//  ParamsCTE.cpp
//  Implementation of the Class ParamsCTE
//  Created on:      05-Dic-2012 16:59:59 PM
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



#include "ParamsCTE.h"



//==============================================================================
/**
 * Constructor method
 */
ParamsCTE::ParamsCTE(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsCTE::~ParamsCTE(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * Generate a CTE map from meanCTE paramater and set this CTE map in the Dataset 
 * to be used in the corresponding processing step.
 */
void ParamsCTE::paramsCTEcreateCTEMap(DataSet &m_DataSet)
{

    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    meanCTE = p_DataSet->datasetGetmeanCTE();
    numLowCTEPixels = p_DataSet->datasetGetnumLowCTEPixels();
    numLowCTELines = p_DataSet->datasetGetnumLowCTELines();
    
    ccdSizeX = p_DataSet->datasetGetccdSizeX();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();

    
    //Generating the CTE map with the CCD size
    cteMap.resize(ccdSizeX, ccdSizeY);
    cteMap = 0.0;
   
    //Assign meanCTE value to each pixel in the map
    cteMap = meanCTE;
    
    //Generate the uniform random distributions randX and randY
    DiscreteUniform<int, ranlib::MersenneTwister, ranlib::independentState> randX(ccdSizeX);
    DiscreteUniform<int, ranlib::MersenneTwister, ranlib::independentState> randY(ccdSizeY);
    randX.seed(DataSet::seedRNG);
    randY.seed(DataSet::seedRNG + 1);
    
    //Number of Low CTE Pixels are assigned with lower CTE value in random positions
    for (int i = 0; i < numLowCTEPixels; i++)
    {
        cteMap(randX.random(), randY.random()) = meanCTE - 50 * (1 - meanCTE); // TODO: The 50 value is HARDCODED. TBC (pma)
    }
    
    //Number of Low CTE lines are assigned with lower CTE value in random rows 
    for (int i = 0; i < numLowCTELines; i++)
    {
        cteMap(randX.random(), Range::all()) = meanCTE - 50 * (1 - meanCTE); // TODO: The 50 value is HARDCODED. TBC (pma)
    }
    
    //Setting the generated CTE map in the Dataset
    m_DataSet.datasetSetCTEMap(cteMap);
    
    
    LogManager::log<< "    CTE map created ";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================
