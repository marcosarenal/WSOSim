///////////////////////////////////////////////////////////
//  PreProcessingPhotometry.h
//  Implementation of the Class PreProcessingPhotometry
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



#ifndef PREPROCESSINGPHOTOMETRY_H_
#define PREPROCESSINGPHOTOMETRY_H_


#include <string>
#include <fstream>
#include <sys/stat.h>

#include "DataSet.h"
#include "ticpp.h"
#include "LogManager.h"
#include "FileUtilities.h"
#include "ParamsPhotometry.h"
#include "DataSet.h"

/**
 * The PreprocessingPhotometry reads all the input files and makes sure that all
 * the photometry parameters to be required by any of the processing steps
 * (included in the ProcessingSteps Subsystem) are included in the DataSet.
 */
class PreProcessingPhotometry
{

public:

	PreProcessingPhotometry();
	PreProcessingPhotometry(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry, std::string photometryParameterFile);        
	virtual ~PreProcessingPhotometry();
	void preProcessingPhotometryReadinput();
    
    
private:

    DataSetPhotometry*  p_datasetPhotometry; //Pointer to the DataSet to retrieve parameters from it.
    DataSet*            p_DataSet;           //Pointer to the DataSet to retrieve parameters from it.

    std::string photometryDirName;                  //Parameter retrieved from DataSet.
    std::string photometryPlotsDir;                   //Output directory for plots.
    std::string outputPath, prefix;               //Parameter retrieved from DataSet.
    std::string photometryMethod;

    int    subPixelsPerPixel;                //Parameter retrieved from DataSet.
    int    psfSubPixels, psfNumPixels;       //Parameter retrieved from DataSet.
    int    psfSize;                          //Parameter retrieved from DataSet.
    int    subFieldSizeX, subFieldSizeY;     //Parameter retrieved from DataSet.

    

    
        

    
};
#endif /* PREPROCESSINGPHOTOMETRY_H_ */ 
