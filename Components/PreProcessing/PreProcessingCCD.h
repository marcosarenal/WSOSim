///////////////////////////////////////////////////////////
//  PreProcessingCCD.h
//  Implementation of the Class PreProcessingCCD
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



#ifndef PREPROCESSINGCCD_H_
#define PREPROCESSINGCCD_H_


#include "ParamsStarfield.h"
#include "LogManager.h"
#include "DataSet.h"
#include "ParamsCCD.h"
#include "ParamsBackground.h"
#include "ParamsCTE.h"
#include "ParamsFlatField.h"
#include "ParamsJitter.h"
#include "ParamsTransit.h"

/**
 * The PreprocessingCCD reads all the input files and makes sure that all the CCD
 * parameters to be required by any of the processing steps (included in the
 * ProcessingSteps Subsystem) are included in the DataSet.
 */
class PreProcessingCCD
{

public:

	PreProcessingCCD(DataSet &m_DataSet);
	virtual ~PreProcessingCCD();
        void preProcessingCCDGetPredefinedCCDOffsets(DataSet *p_DataSet);


private:
    
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
        
	int     ccdSizeX, ccdSizeY;                     //Parameter retrieved from DataSet.
	double  pixelSize;                              //Parameter retrieved from DataSet.
	double  originOffsetXmm, originOffsetYmm;       //Parameter retrieved from DataSet.
	double  ccdOrientation;                         //Parameter retrieved from DataSet.
        
	std::string  convolutionMethod;                      //Parameter retrieved from DataSet.
	std::string  ccdPredefinedPosition;                  //Parameter retrieved from DataSet.
	std::string  parameterFile;                          //Parameter retrieved from DataSet.
	bool    performExoTransit;                      //Parameter retrieved from DataSet.
        
        

};
#endif /* PREPROCESSINGCCD_H_ */
