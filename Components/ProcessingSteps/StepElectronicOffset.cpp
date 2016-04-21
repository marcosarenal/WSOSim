///////////////////////////////////////////////////////////
//  StepElectronicOffset.cpp
//  Implementation of the Class StepElectronicOffset
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



#include "StepElectronicOffset.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepElectronicOffset::StepElectronicOffset(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepElectronicOffset::~StepElectronicOffset(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies electronic offset to the pixelMap, smearingMap and biasRegisterMap.
 */
void StepElectronicOffset::StepElectronicOffsetapplication(DataSet &m_DataSet)
{
    
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	

    //Retrieving parameters from DataSet
    electronicOffset = p_DataSet->datasetGetelectronicOffset();     


    //Retrieving the BIAS register map from DataSet        
    //Initialize smearingMap
    biasRegisterMap.resize(m_DataSet.datasetGetbiasRegisterMap().extent(0), m_DataSet.datasetGetbiasRegisterMap().extent(1));
    biasRegisterMap = 0.0;
    biasRegisterMap = m_DataSet.datasetGetbiasRegisterMap();

         
    //Initialize smearingMap
    smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
    smearingMap = 0.0;
    smearingMap = m_DataSet.datasetGetsmearingMap();

    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();


    pixelMap += electronicOffset;

    biasRegisterMap += electronicOffset;

    smearingMap += electronicOffset;
                
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);
    
    //Setting the new calculated biasRegisterMap into the DataSet
    m_DataSet.datasetSetbiasRegisterMap(biasRegisterMap);
    
    //Setting the new calculated smearingMap into the DataSet
    m_DataSet.datasetSetsmearingMap(smearingMap);



    LogManager::log <<"    Successfully applied Electronic Offset.";
    GlobalVariables::logManager.LogManagerShowLog();     

}
//=========================================================================