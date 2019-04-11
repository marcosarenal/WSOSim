///////////////////////////////////////////////////////////
//  StepGain.h
//  Implementation of the Class StepGain
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



#ifndef STEPGAIN_H_
#define STEPGAIN_H_



#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"



/**
 * Class for applying a given gain to the pixels of a CCD/CMOS. Only the normal
 * pixels, prescan and overscan of a CCD object are affected.
 */
class StepGain
{

public:
	StepGain();
	virtual ~StepGain();
        void StepGainapplication(DataSet &m_DataSet);

        
private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
        double gain;                                                    //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap, smearingMap, biasRegisterMap;        //Map array retrieved from DataSet.

};
#endif /* STEPGAIN_H_ */