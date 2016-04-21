///////////////////////////////////////////////////////////
//  StepSaturation.h
//  Implementation of the Class StepSaturation
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



#ifndef STEPSATURATION_H_
#define STEPSATURATION_H_


#include "DataSet.h"
#include "blitz/array.h"
#include "FileUtilities.h"

using namespace blitz;



/**
 * Compute the effects of saturation on the detector. If a pixel is over-exposed,
 * the extra electrons are distributed evenly on the neighboring pixels in readout
 * direction. At the edges of the detector (here actually the sub-image) the flux
 * is lost. Digital saturation results in a cutoff of pixel flux that is higher
 * than the A/D converter can handle (e.g. on a 16-bit converter a flux value of
 * 70000 is changed to 65535 (=2^16-1). The saturation values are taken from the
 * CCD object.
 */
class StepSaturation
{

public:
	StepSaturation();
	virtual ~StepSaturation();
        void StepSaturationfullWellapplication (DataSet &m_DataSet);
        void StepSaturationdigitalapplication (DataSet &m_DataSet);

        
        private:

        DataSet*    p_DataSet;               //Pointer to the DataSet to retrieve parameters from it.
    
        int    fullWellSat;                  //Parameter retrieved from DataSet.
        int    digitalSat;                   //Parameter retrieved from DataSet.

        Array<float, 2>  pixelMap;           //Map array retrieved from DataSet.
        Array<float, 2>  smearingMap;        //Map array retrieved from DataSet.
        Array<float, 2>  biasRegisterMap;    //Map array retrieved from DataSet.

};
#endif /* STEPSATURATION_H_ */