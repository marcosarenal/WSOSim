///////////////////////////////////////////////////////////
//  StepRebin.h
//  Implementation of the Class StepRebin
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



#ifndef STEPREBIN_H_
#define STEPREBIN_H_

#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"


using namespace blitz;


/**
 * Class for applying the rebinning, in the case that the subpixelMap has been used 
 * instead of the pixelMap.
 */
class StepRebin
{

public:
	StepRebin();
	virtual ~StepRebin();
        void StepRebinapplication(DataSet &m_DataSet);
        
private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
        int    subPixelsPerPixel, subFieldSizeX, subFieldSizeY;         //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap, subPixelMap;                         //Map array retrieved from DataSet.
};
#endif /* STEPREBIN_H_ */