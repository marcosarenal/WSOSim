///////////////////////////////////////////////////////////
//  PreProcessingPSF.cpp
//  Implementation of the Class PreProcessingPSF
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




#include "PreProcessingPSF.h"



//==============================================================================
/**
 * Constructor method
 */
PreProcessingPSF::PreProcessingPSF(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
PreProcessingPSF::~PreProcessingPSF(){}
//==============================================================================



//==============================================================================
/**
 * The Preprocessing PSF triggers the ParamsPSF module in charge of calculate the 
 * PSF matrix and Set it in the Dataset. It is done separately from the PreprocessingCCD
 * module because some other parameters calculated in ParamsPSF and set in the Dataset
 * needs to be calculated in advance of the PreprocessingCCD.
 */
PreProcessingPSF::PreProcessingPSF(DataSet &m_DataSet)
{
 
    //Calculate the PSF parameters to be set in the DataSet.
    m_ParamsPSF.paramsPSFCalculation(m_DataSet);

}
//==============================================================================


