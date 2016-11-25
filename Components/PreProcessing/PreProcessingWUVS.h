///////////////////////////////////////////////////////////
//  PreProcessingWUVS.h
//  Implementation of the Class PreProcessingWUVS
//  Created on:      11-June-2015 1:59:59 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
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



#ifndef PREPROCESSINGWUVS_H_
#define PREPROCESSINGWUVS_H_


#include "ParamsStarfield.h"
#include "LogManager.h"
#include "DataSet.h"
#include "ParamsWUVS.h"
#include "ParamsBackground.h"
#include "ParamsCTE.h"
#include "ParamsFlatField.h"
#include "ParamsJitter.h"

/**
 * The PreprocessingCCD reads all the input files and makes sure that all the CCD
 * parameters to be required by any of the processing steps (included in the
 * ProcessingSteps Subsystem) are included in the DataSet.
 */
class PreProcessingWUVS
{

public:

	PreProcessingWUVS(DataSet &m_DataSet);
	virtual ~PreProcessingWUVS();
        void preProcessingCCDGetPredefinedCCDOffsets(DataSet *p_DataSet);


private:
    
        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
        
	int     ccdSizeX, ccdSizeY;                     //Parameter retrieved from DataSet.
	double  pixelSize;                              //Parameter retrieved from DataSet.
	double  originOffsetXmm, originOffsetYmm;       //Parameter retrieved from DataSet.
	double  ccdOrientation;                         //Parameter retrieved from DataSet.
        
	std::string  convolutionMethod;                      //Parameter retrieved from DataSet.
	std::string  ccdPredefinedPosition;                  //Parameter retrieved from DataSet.
	std::string parameterFile;                           //Parameter retrieved from DataSet.
            
        Array<float, 2> inputImage;
        Array<int, 2> inputIntImage;
        std::string  outputPath, prefix;                               //Parameter retrieved from DataSet.

        

};
#endif /* PREPROCESSINGWUVS_H_ */
