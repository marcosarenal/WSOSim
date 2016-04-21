///////////////////////////////////////////////////////////
//  ParamsCTE.h
//  Implementation of the Class ParamsCTE
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



#ifndef PARAMSCTE_H_
#define PARAMSCTE_H_

#include "blitz/array.h"
#include "random/discrete-uniform.h"
#include "DataSet.h"
//#include "GlobalVariables.h"

using namespace ranlib;




/**
 * This class creates a CTE map that will be used in the StepCTE component.
 */
class ParamsCTE
{
    
public:
	ParamsCTE();
	virtual ~ParamsCTE();    
    void paramsCTEcreateCTEMap(DataSet &m_DataSet);
    
private:
    
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    double meanCTE;                                 //Parameter retrieved from DataSet.
    int    numLowCTEPixels, numLowCTELines;         //Parameter retrieved from DataSet.
    int    ccdSizeX, ccdSizeY;                      //Parameter retrieved from DataSet.
    
    Array<double, 2> cteMap;                        //Blitz map set into the DataSet.
    
};
#endif /* PARAMSCTE_H_ */

