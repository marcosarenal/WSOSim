///////////////////////////////////////////////////////////
//  ParamsFlatField.h
//  Implementation of the Class ParamsFlatField
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



#ifndef PARAMSFLATFIELD_H_
#define PARAMSFLATFIELD_H_

#include "blitz/array.h"
#include "random/discrete-uniform.h"
#include "random/normal.h"
#include "DataSet.h"
//#include "GlobalVariables.h"
#include "FFT.h"


using namespace blitz;
using namespace ranlib;




/**
 * This class creates a FlatField map that will be used in the StepFlatField component.
 */
class ParamsFlatField
{
    
public:
	ParamsFlatField();
	virtual ~ParamsFlatField();
    void paramsFlatFieldcreateFFMap(DataSet &m_DataSet);
    void paramsFlatFieldcreateSubPixelFFMap();
    void paramsFlatFieldcreatePixelFFMap();
    
    
    
private:
    
    //Pointer to DataSet to retrieve parameters from it.
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
    
    
	double flatfieldPixelNoise, flatfieldWhiteNoise;        //Parameter retrieved from DataSet.
	double flatfieldIntraPixelWidth;                        //Parameter retrieved from DataSet.
    int    pixelSize, subPixelsPerPixel;                    //Parameter retrieved from DataSet.
    int    edge;
    int    subFieldSizeX, subFieldSizeY;                    //Parameter retrieved from DataSet.
    
    Array<float, 2>  flatfieldMap, subflatfieldMap;         //Blitz map set into the DataSet.
    Array<double, 1> intraPixelSensitivity;
    
    
};
#endif /* PARAMSFLATFIELD_H_ */

