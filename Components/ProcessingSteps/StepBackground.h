///////////////////////////////////////////////////////////
//  StepBackground.h
//  Implementation of the Class StepBackground
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



#ifndef STEPBACKGROUND_H_
#define STEPBACKGROUND_H_


#include "LogManager.h"
#include "DataSet.h"

/**
 * Computes and sets the sky background value on the CCD/CMOS for a certain
 * location in the sky.
 */
class StepBackground
{

public:
	StepBackground();
	virtual ~StepBackground();
        void StepBackgroundapplication(DataSet &m_DataSet);

private:

        DataSet     m_DataSet;                         //m_DataSet object declaration. 
        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
        

        int    subFieldSizeX, subPixelsPerPixel;        //Parameter retrieved from DataSet.
        double background, exposureTime;                //Parameter retrieved from DataSet.

        Array<float, 2>  subPixelMap;                   //Map array retrieved from DataSet.
 
};
#endif /* STEPBACKGROUND_H_ */
