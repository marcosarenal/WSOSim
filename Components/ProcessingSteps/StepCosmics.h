///////////////////////////////////////////////////////////
//  StepCosmics.h
//  Implementation of the Class StepCosmics
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



#ifndef STEPCOSMICS_H_
#define STEPCOSMICS_H_

#include "LogManager.h"
#include "DataSet.h"
#include "Statistics.h"
#include "blitz/array.h"
#include "random/discrete-uniform.h"
#include "Constants.h"

/**
 * Class for producing stochastically distributed (in time and space) cosmic hits
 * on the CCD/CMOS. Mean number of cosmic hits in units of events/cm^2/minute. The
 * number of occurrences on the detector sub-field is randomly Poisson distributed
 * and depends on the parameter hitRate, the size of the sub-field, and the
 * exposure time. The intensity and shape of the simulated cosmic hits and are
 * normally distributed random values. Cosmics are simulated through a 2-
 * dimensional general (non-circular) Gaussian function. The flux of the central
 * normal pixel of the Gaussian is indicated through the parameters
 * CosmicSaturation in units of the full well saturation. The width of the
 * Gaussian in normal pixels is defined by the parameter CosmicWidth. The cosmics
 * only affect the simulated detector sub-field. Neither, the simulated trailing-
 * row, flatfield, or overscan are affected.
 */
class StepCosmics
{

public:
	StepCosmics();
	virtual ~StepCosmics();
        void StepCosmicsapplication(DataSet &m_DataSet);
        void StepCosmicsadd2DGauss(int centerX, int centerY, double intensity, double sigmax, double sigmay, double theta);

        
private:

        DataSet*    p_DataSet;                       //Pointer to the DataSet to retrieve parameters from it.
        
        int     subPixelsPerPixel;                   //Parameter retrieved from DataSet.
        int     fullWellSat;                         //Parameter retrieved from DataSet.
        int     subFieldSizeX, subFieldSizeY;        //Parameter retrieved from DataSet.
        double  exposureTime;                        //Parameter retrieved from DataSet.
        double  pixelSize;                           //Parameter retrieved from DataSet.
        double  cosmicHitRate;                       //Parameter retrieved from DataSet.
        double  cosmicsWidth;                        //Parameter retrieved from DataSet.
        double  cosmicsLength;                       //Parameter retrieved from DataSet.

        double  cosmicsSatFactor;                    //Parameter retrieved from DataSet.

        double  ccdSubFieldSizeSqcm; 
        double  meanEvents; 


        Array<float, 2>  subPixelMap;                //Map array retrieved from DataSet.


};
#endif /* STEPCOSMICS_H_ */