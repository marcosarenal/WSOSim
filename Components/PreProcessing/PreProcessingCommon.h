///////////////////////////////////////////////////////////
//  PreProcessingCommon.h
//  Implementation of the Class PreProcessingCommon
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



#ifndef PREPROCESSINGCOMMON_H_
#define PREPROCESSINGCOMMON_H_

#include <string>
#include <fstream>
#include <sys/stat.h>

#include "DataSet.h"
#include "ticpp.h"
#include "LogManager.h"
#include "FileUtilities.h"



using namespace std;
using namespace ticpp;

/**
 * This class performs the duties of the PreProcessing relating to the common parameters of any processing mode.
 * This includes calling the datasetReadParameterFile(parameterFile) method to read the input parameters file 
 * and set all these parameters in the DataSet.
 */
class PreProcessingCommon
{

public:
        
	PreProcessingCommon();
	virtual ~PreProcessingCommon();
        
	void PreProcessingCommonCCD(DataSet &m_DataSet, std::string parameterFile);
        void PreProcessingCommonPhotometry(DataSetPhotometry &datasetPhotometry, std::string photometryParameterFile);

	     
        

private:

        DataSet*    p_DataSet;                  //Pointer to the DataSet to retrieve parameters from it.

     //   std::string parameterFile;                   //Parameter retrieved from DataSet.
        std::string outputDir;                       //Parameter retrieved from DataSet.
        std::string outputPath, prefix;              //Parameter retrieved from DataSet.
        
};
#endif /* PREPROCESSINGCOMMON_H_ */
