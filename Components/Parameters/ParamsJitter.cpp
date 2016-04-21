///////////////////////////////////////////////////////////
//  ParamsJitter.cpp
//  Implementation of the Class ParamsJitter
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



#include "ParamsJitter.h"



//==============================================================================
/**
 * Constructor method
 */
ParamsJitter::ParamsJitter(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsJitter::~ParamsJitter(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * Sets the Jitter time, yaw, pitch and roll input parameters in the Dataset to 
 * be used in the corresponding processing step.
 */
void ParamsJitter::ParamsJittersetParameters(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    jitterFile = p_DataSet->datasetGetjitterFile();
    jitterMultFactor = p_DataSet->datasetGetjitterMultFactor();
	
    //Defining Blitz arrays containing yaw, pitch, roll and jitterTimes
    Array<double, 1> yaw;
    Array<double, 1> pitch;
    Array<double, 1> roll;
    Array<double, 1> jitterTimes;
    
    
    //Generating the jitter input array
    FileUtilities::readCol(jitterFile, 0, jitterTimes);
    FileUtilities::readCol(jitterFile, 1, yaw);
    FileUtilities::readCol(jitterFile, 2, pitch);
    FileUtilities::readCol(jitterFile, 3, roll);

    yaw = yaw * jitterMultFactor;
    pitch = pitch * jitterMultFactor;
    roll = roll * jitterMultFactor;
    
    jitterInputParams.resize(yaw.rows(), 4);
    jitterInputParams = 0.0;
    jitterInputParams(Range::all(), 0) = jitterTimes(Range::all());
    jitterInputParams(Range::all(), 1) = yaw(Range::all());
    jitterInputParams(Range::all(), 2) = pitch(Range::all());
    jitterInputParams(Range::all(), 3) = roll(Range::all());
    
    
    //Setting the generated jitterInputParams map in the Dataset
    m_DataSet.datasetSetJitterInputParams(jitterInputParams);
    
    
    //Free memory
    jitterTimes.free();
    yaw.free();
    pitch.free();
    roll.free();
    
    
    LogManager::log<< "    Jitter parameters set ";
    GlobalVariables::logManager.LogManagerShowLog();

}
//==============================================================================
