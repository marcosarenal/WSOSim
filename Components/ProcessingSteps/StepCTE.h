///////////////////////////////////////////////////////////
//  StepCTE.h
//  Implementation of the Class StepCTE
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



#ifndef STEPCTE_H_
#define STEPCTE_H_

#include "blitz/array.h"
#include "FileUtilities.h"
#include "LogManager.h"
#include "DataSet.h"


/**
 * Simulate charge transfer efficiency degradation on a CCD.
 */
class StepCTE
{

public:
	StepCTE();
	virtual ~StepCTE();
        void StepCTEapplication(DataSet &m_DataSet);


private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.

        int    subFieldZeroX, subFieldZeroY;            //Parameter retrieved from DataSet.
        int    subFieldSizeX, subFieldSizeY;            //Parameter retrieved from DataSet.
        double meanCTE;                                 //Parameter retrieved from DataSet.
        int    numLowCTEPixels, numLowCTELines;         //Parameter retrieved from DataSet.
        int    ccdSizeY;                                //Parameter retrieved from DataSet.
        double exposureTime;                            //Parameter retrieved from DataSet.
        double readOutTime;                             //Parameter retrieved from DataSet.

        double rowReadOutTime;                          //Time from the beginning of the CCD read out to reading a determined row 

        Array<double, 2> cteMap;                        //Map array retrieved from DataSet.
        Array<float, 2> pixelMap;                       //Map array retrieved from DataSet.



};
#endif /* STEPCTE_H_ */
