///////////////////////////////////////////////////////////
//  PreProcessingCMOS.h
//  Implementation of the Class PreProcessingCMOS
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



#ifndef PREPROCESSINGCMOS_H_
#define PREPROCESSINGCMOS_H_

#include "ParamsCMOS.h"
#include "DataSet.h"

/**
 * The PreprocessingCMOS reads all the input files and makes sure that all CMOS
 * parameters to be required by any of the processing steps (included in the
 * ProcessingSteps Subsystem) are included in the DataSet.
 */
class PreProcessingCMOS
{

public:
	PreProcessingCMOS();
	virtual ~PreProcessingCMOS();
	void preProcessingCMOSParams();
	void preProcessingCMOSReadinput();
private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
};
#endif /* PREPROCESSINGCMOS_H_ */
