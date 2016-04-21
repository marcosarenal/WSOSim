///////////////////////////////////////////////////////////
//  PostProcessing.h
//  Implementation of the Class PostProcessing
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



#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

#include "DataSet.h"
#include "FileUtilities.h"
#include "CCfits/CCfits"
#include "blitz/array.h"

using namespace CCfits;

/**
 * Convenience functions to write files.
 */
class PostProcessing
{
    
public:
    
	PostProcessing();
	virtual ~PostProcessing();
    PostProcessing(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry);
    void postProcessingGenerateNoisePlots();

    
    //	void postProcessingWriteCCDParams();
    //	void postProcessingWriteCMOSParams();
	void postProcessingWriteXML();
    
private:
    
    DataSet* p_DataSet;                                     //Pointer to the DataSet to retrieve parameters from it.
    DataSetPhotometry* p_datasetPhotometry;                 //Pointer to the DataSet to retrieve parameters from it.
    string fileName;                                        //Parameter retrieved from DataSet.
    string outputPath;                                      //Parameter retrieved from DataSet.
    string prefix;                                          //Parameter retrieved from DataSet.
    string outputDir;                                       //Parameter retrieved from DataSet.
    string mapName;                                         //Name of the Blitz array to be written to FITS file.
    string photometryPlotsDir;                              //Parameter retrieved from DataSetPhotometry.
    string photometryDirName;                               //Parameter retrieved from DataSetPhotometry.
    double background;                                      //Parameter retrieved from DataSetPhotometry.
    int    photometryNumTelescopes;                         //Parameter retrieved from DataSetPhotometry.

    double exposureTime;                                    //Parameter retrieved from DataSet.
    double flatfieldIntraPixelWidth;                        //Parameter retrieved from DataSet.
    
    Array<float, 2> pixelMap, psfMap, subPixelMap;          //Array retrieved from DataSet.
    Array<float, 2>  flatfieldMap, smearingMap;             //Array retrieved from DataSet.
    Array<double, 2> cteMap;                                //This is a Blitz 2-D array with the CTE map.
    
    
};
#endif /* POSTPROCESSING_H_ */
