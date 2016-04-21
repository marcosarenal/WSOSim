///////////////////////////////////////////////////////////
//  StepPhotonNoise.h
//  Implementation of the Class StepPhotonNoise
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



#ifndef STEPPHOTONNOISE_H_
#define STEPPHOTONNOISE_H_

#include "DataSet.h"
#include "blitz/array.h"
#include "Statistics.h"
#include "FileUtilities.h"


using namespace blitz;



/**
 * Adds photon (shot) noise to the image of a detector object.
 */
class StepPhotonNoise
{

public:
	StepPhotonNoise();
	virtual ~StepPhotonNoise();
        void StepPhotonNoiseapplication(DataSet &m_DataSet);
        void StepPhotonNoiseWUVSapplication(DataSet &m_DataSet);
        
        

private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.

        int    numSmearingOverscanRows;                //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap;                     //Map array retrieved from DataSet.
        Array<float, 2>  smearingMap;                  //Map array retrieved from DataSet.

};


#endif /* STEPPHOTONNOISE_H_ */