///////////////////////////////////////////////////////////
//  StepReadOutNoise.h
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



#ifndef STEPREADOUTNOISE_H_
#define STEPREADOUTNOISE_H_


#include "DataSet.h"
#include "blitz/array.h"
#include "Statistics.h"
#include "random/normal.h"
#include "random/uniform.h"


using namespace ranlib;
using namespace std;
using namespace blitz;


/**
 * Compute the effects of readout noise for a CCD/CMOS image on pixel level. The
 * mean readout noise is a parameter of CCD/CMOS. We add a normally distributed
 * value around the mean readout noise value to the image.
 */
class StepReadOutNoise
{

public:
	StepReadOutNoise();
	virtual ~StepReadOutNoise();
        void StepReadOutNoiseapplication(DataSet &m_DataSet);

        
private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
	double readOutNoise;                           //Parameter retrieved from DataSet.
        int    numPrescanRows;                         //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap;                     //Map array retrieved from DataSet.
        Array<float, 2>  biasRegisterMap;              //Map array retrieved from DataSet.

        

};
#endif /* STEPREADOUTNOISE_H_ */