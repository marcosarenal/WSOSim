///////////////////////////////////////////////////////////
//  StepBackground.cpp
//  Implementation of the Class StepBackground
//  Created on:      23-Oct-2012 2:00:00 PM
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



#include "StepBackground.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepBackground::StepBackground(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepBackground::~StepBackground(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 *This function adds the Background effect to the image.
 */
void StepBackground::StepBackgroundapplication(DataSet &m_DataSet)
{

    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    //Retrieving parameters from DataSet
    exposureTime = p_DataSet->datasetGetexposureTime();

    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    background = p_DataSet->datasetGetbackground(); //Given in e/(s*pixel)
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();

    //Initializing and retrieving the subpixel map from DataSet
    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
    subPixelMap = 0.0;
    subPixelMap = m_DataSet.datasetGetsubPixelMap();


    //Setting the background in the subPixelMap
    subPixelMap += (background * exposureTime / (subPixelsPerPixel * subPixelsPerPixel));



    //Setting the new calculated subPixelMap into the DataSet
    m_DataSet.datasetSetsubPixelMap(subPixelMap);  


    LogManager::log <<"    Successfully added Background effect.";
    GlobalVariables::logManager.LogManagerShowLog();     



}
//==============================================================================