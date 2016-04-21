///////////////////////////////////////////////////////////
//  ParamsJitter.h
//  Implementation of the Class ParamsJitter
//  Created on:      23-Oct-2012 1:59:59 PM
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



#ifndef PARAMSJITTER_H_
#define PARAMSJITTER_H_

#include "blitz/array.h"
#include "DataSet.h"
#include "FileUtilities.h"

using namespace blitz;




/**
 * This class sets the Jitter time, yaw, pitch and roll input parameters in the Dataset to 
 * be used in the StepJitter component.
 */
class ParamsJitter
{

public:
	ParamsJitter();
	virtual ~ParamsJitter();
        void ParamsJittersetParameters(DataSet &m_DataSet);

private:
    
        DataSet* p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.



        string  jitterFile;                                     //Parameter retrieved from DataSet.
        double  jitterMultFactor;                               //Parameter retrieved from DataSet.
	    Array<double, 2> jitterInputParams;                     //Generated Blitz map to be set into the DataSet.     

};
#endif /* PARAMSJITTER_H_ */

