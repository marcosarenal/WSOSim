///////////////////////////////////////////////////////////
//  StepFlatField.h
//  Implementation of the Class StepFlatField
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



#ifndef STEPFLATFIELD_H_
#define STEPFLATFIELD_H_

#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"


using namespace blitz;


/**
 * Class for computing the Pixel Non-Uniform Response or flatfield map. The
 * flatfield ( pixel non-uniform response) is computed by considering a spatial
 * 1/f-response of the sensitivity.
 */
class StepFlatField
{

public:
	StepFlatField();
	virtual ~StepFlatField();
        void StepFlatFieldapplicationSubpixel(DataSet &m_DataSet);
        void StepFlatFieldapplicationPixel(DataSet &m_DataSet);

private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.

        Array<float, 2>  pixelMap, subPixelMap;                 //Map array retrieved from DataSet.
        Array<float, 2>  flatfieldMap, subflatfieldMap;         //Map array retrieved from DataSet.
};
#endif /* STEPFLATFIELD_H_ */