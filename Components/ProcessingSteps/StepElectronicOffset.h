///////////////////////////////////////////////////////////
//  StepElectronicOffset.h
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



#ifndef STEPELECTRONICOFFSET_H_
#define STEPELECTRONICOFFSET_H_


#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"

/**
 * Add a fixed electronic offset value to the CCD/CMOS. Affects the normal pixel
 * map (NOT the sub-pixel map), the pre-scan and the over-scan.
 */
class StepElectronicOffset
{

public:
	StepElectronicOffset();
	virtual ~StepElectronicOffset();
        void StepElectronicOffsetapplication(DataSet &m_DataSet);

private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
        double electronicOffset;                                //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap, smearingMap, biasRegisterMap;        //Map array retrieved from DataSet.

};
#endif /* STEPELECTRONICOFFSET_H_ */