///////////////////////////////////////////////////////////
//  StepReadOutNoise.cpp
//  Implementation of the Class StepReadOutNoise
//  Created on:      23-Oct-2012 2:00:02 PM
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



#include "StepReadOutNoise.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepReadOutNoise::StepReadOutNoise(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepReadOutNoise::~StepReadOutNoise(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the Read Out Noise effect to the pixelMap.
 */
void StepReadOutNoise::StepReadOutNoiseapplication(DataSet &m_DataSet)
{
    
        //Pointing to the DataSet
        p_DataSet = &m_DataSet;       	
        
        //Retrieving parameters from DataSet
        readOutNoise = p_DataSet->datasetGetreadOutNoise();
        numPrescanRows = p_DataSet->datasetGetnumPrescanRows();
                
                
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        pixelMap = m_DataSet.datasetGetpixelMap();


        //Create a BIAS register Map and initialize it.
        biasRegisterMap.resize(pixelMap.rows(), numPrescanRows);
        biasRegisterMap = 0.0;

        //Generate a Normal distribution Map called randomMap with mean value = readOutNoise and standard Deviation = sqrt(readOutNoise)
        Normal<double, ranlib::MersenneTwister, ranlib::independentState> randomMap(readOutNoise, sqrt(readOutNoise));
        randomMap.seed(DataSet::seedRNG);
        
        //For each of the pixels in the biasRegisterMap, added a random value taken from the randomMap
        for (int i = 0; i < pixelMap.rows(); i++)
        {
            for (int j = 0; j < numPrescanRows; j++)
            {
                biasRegisterMap(i, j) += randomMap.random();
            }
        }
        
        //For each of the pixels in the pixelMap, added a random value taken from the randomMap
        for (int i = 0; i < pixelMap.rows(); i++)
        {
            for (int j = 0; j < pixelMap.cols(); j++)
            {
                pixelMap(i, j) += randomMap.random();
            }
        }


                    
        //Setting the new calculated pixelMap into the DataSet
        m_DataSet.datasetSetpixelMap(pixelMap);

        //Setting the new calculated biasRegisterMap into the DataSet
        m_DataSet.datasetSetbiasRegisterMap(biasRegisterMap);

        
                        
        LogManager::log <<"    Successfully added Read Out Noise effect.";
        GlobalVariables::logManager.LogManagerShowLog();     


}
//==============================================================================

