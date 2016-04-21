///////////////////////////////////////////////////////////
//  StepCTE.cpp
//  Implementation of the Class StepCTE
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



#include "StepCTE.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepCTE::StepCTE(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepCTE::~StepCTE(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 *This function adds the charge transfer efficiency effect to the image.
 */
void StepCTE::StepCTEapplication(DataSet &m_DataSet)
{

    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Initialize cteMap
    cteMap.resize(m_DataSet.datasetGetCTEMap().extent(0), m_DataSet.datasetGetCTEMap().extent(1));
    cteMap = 0.0;
    //Retrieving the cteMap from DataSet
    cteMap = m_DataSet.datasetGetCTEMap();


    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();
    
    //Retrieving parameters from DataSet
    subFieldZeroX  = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY  = p_DataSet->datasetGetsubFieldZeroY();
    subFieldSizeX  = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY  = p_DataSet->datasetGetsubFieldSizeY();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();
    
    exposureTime = p_DataSet->datasetGetexposureTime();
    readOutTime = p_DataSet->datasetGetreadOutTime();

        
        
    //The complete (true size) CCD is read out
    Array<float, 2> shiftMap(subFieldSizeX, subFieldZeroY + subFieldSizeY);
    shiftMap = 0.;
    //copy the sub image into the shiftMap
    shiftMap(Range::all(), Range(subFieldZeroY, subFieldZeroY + subFieldSizeY - 1)) = pixelMap;
    Array<float, 2> shiftMapTemp(shiftMap.shape());


    Array<float, 2> newMap(subFieldSizeX, subFieldZeroY + subFieldSizeY);
    Array<double, 1> readoutStrip(subFieldSizeX);
    readoutStrip = 0;


    //we assume that the serial register has a CTE of 1, unlike the CCD that have a cteMap
    
    //Read out the CCD in columns (in the Y axis), row by row
    for (int rowIter = 0; rowIter < subFieldZeroY + subFieldSizeY; rowIter++) //for each row in the CCD
    {
        //Calculate the time required to for reading a certain row (rowIter) as the total read out time divided by the number of rows(subFieldSizeY).  
        rowReadOutTime = rowIter * readOutTime/ccdSizeY;

        //Copy the bottom row to the readout register multiplied by its corresponding row in the cteMap. 
        //The read out register encounters the flux of the previously read rows as a function of the time required to read that row.
        readoutStrip += cteMap(Range(subFieldZeroX, subFieldZeroX + subFieldSizeX - 1), 0) * (shiftMap(Range::all(), 0)) * rowReadOutTime/ exposureTime;


        //For each pixel in the current row
        for (int pixelIter = 0; (unsigned)pixelIter < readoutStrip.size(); pixelIter++) 
        {
            //Add the pixel value to the map containing the charge transfer contamination 
            newMap(pixelIter, rowIter) = readoutStrip(pixelIter) ;
        }

        //Shift the electrons of CCD towards readout register: Shift the map one row down
        for (int x = 0; x < subFieldSizeX; x++)
        {
            for (int y = 0; y < subFieldZeroY + subFieldSizeY - 1; y++) 
            {
                //For each pixel in the shiftedMap
                //its value is obtained by removing the previous value (1-shiftMap(x, y)) and adding the value of the pixel above itself (shiftMap(x, y + 1)) 
                //both multiplied by the cteMap
                shiftMapTemp(x, y) = (1 - cteMap(x + subFieldZeroX, y)) * shiftMap(x, y) + cteMap(x + subFieldZeroX, y + 1) * shiftMap(x, y + 1);

            }
        }

        // New position of the  shiftMap
        shiftMap = shiftMapTemp;
    }

    //Add the charge transfer contamination to the pixelMap
    pixelMap += newMap(Range::all(), Range(subFieldZeroY, subFieldZeroY + subFieldSizeY - 1));
    
                
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);



    //Free temporary maps
    shiftMap.free();
    newMap.free();
    readoutStrip.free();
    shiftMapTemp.free();
        
        
    LogManager::log <<"    Successfully added Charge Transfer contamination effect (CTE).";
    GlobalVariables::logManager.LogManagerShowLog();    
}

//==============================================================================
