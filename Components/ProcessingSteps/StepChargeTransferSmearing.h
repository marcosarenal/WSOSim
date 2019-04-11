///////////////////////////////////////////////////////////
//  StepChargeTransferSmearing.h
//  Implementation of the Class StepChargeTransferSmearing
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



#ifndef STEPCHARGETRANSFERSMEARING_H_
#define STEPCHARGETRANSFERSMEARING_H_

#include "DataSet.h"
#include "FFT.h"
#include "LogManager.h"
#include "FileUtilities.h"



/**
 * Simulate the flux smearing effect caused by the transfer of the charge while
 * the shutter is open.
 */
class StepChargeTransferSmearing
{

public:
	StepChargeTransferSmearing();
	virtual ~StepChargeTransferSmearing();
        void StepChargeTransferSmearingapplication(DataSet &m_DataSet);
        void StepChargeTransferSmearinggenerateSmearingMap(Array<float, 1> &conv);
        void StepChargeTransferSmearingconvolve(Array<float, 1> data, Array<float, 1> mask, Array<float, 1> &conv);
        void StepChargeTransferSmearingremapArrays(Array<float, 1> &data, Array<float, 1> &mask);


private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
        
        int    subFieldSizeX, subFieldSizeY, subPixelsPerPixel;         //Parameter retrieved from DataSet.
        double readOutTime;                                             //Parameter retrieved from DataSet.
        double exposureTime;                                            //Parameter retrieved from DataSet.
        int    numSmearingOverscanRows;                                 //Parameter retrieved from DataSet.
        int    subFieldZeroX;                                           //Parameter retrieved from DataSet.
        double background;                                              //Parameter retrieved from DataSet.
        double fluxm0;                                                  //Parameter retrieved from DataSet.
        double areaTelescope, transEff, quantEff;                       //Parameter retrieved from DataSet.
        int    ccdSizeY;                                                //Parameter retrieved from DataSet.
        int    edgePixels;                                              //Parameter retrieved from DataSet.


 
	    Array<float, 2> smearingMap;            //Map array to be set into the DataSet.
        Array<float, 2> psfMap;                 //Map array retrieved from DataSet.
        Array<float, 2> pixelMap;               //Map array retrieved from DataSet.
        Array<float, 2> subPixelMap;            //Map array retrieved from DataSet.
        Array<float, 2> starListOnCCD;          //Map array retrieved from DataSet.

};
#endif /* STEPCHARGETRANSFERSMEARING_H_ */