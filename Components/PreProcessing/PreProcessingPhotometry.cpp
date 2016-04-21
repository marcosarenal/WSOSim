///////////////////////////////////////////////////////////
//  PreProcessingPhotometry.cpp
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



#include "PreProcessingPhotometry.h"




//==============================================================================
/**
 * Constructor method
 */
PreProcessingPhotometry::PreProcessingPhotometry(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
PreProcessingPhotometry::~PreProcessingPhotometry(){}
//==============================================================================



//==============================================================================
/**
 * The PreProcessingCommon is in charge of reading the input file and creating the output
 * directory using the path and name set in this input file; 
 * @param m_DataSet DataSet object designed to store all the input parameters (among others)
 *        included into the xml parameters input file.
 * @param photometryParameterFile Name of the xml photometry parameters input file.
 */
PreProcessingPhotometry::PreProcessingPhotometry(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry, string photometryParameterFile)
{    
    //Point to the DataSet
    p_datasetPhotometry = &datasetPhotometry;
    p_DataSet = &m_DataSet;
    
            
    //Read parameters file.
    p_datasetPhotometry->datasetPhotometryReadParameterFile(photometryParameterFile);
    
    //Check subfield size input parameter for Weighted Mask photometry 
    photometryMethod = p_datasetPhotometry->datasetGetphotometryMethod();
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels();
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();

    //If PSFsubpixels < subfield size
    if (photometryMethod == "WM" && (psfSubPixels < subFieldSizeX || psfSubPixels < subFieldSizeY) )
    {
        cerr << "\nError (PreProcessingPhotometry::PreProcessingPhotometry()):  You must provide a higher PSFSubPixels value (being PSFSubPixels/SubPixels an integer value) or a smaller subFieldSize because" << endl;
        cerr << "psfSubPixels ("<< psfSubPixels <<")  < subFieldSizeX ("<< subFieldSizeX <<") or psfSubPixels ("<< psfSubPixels <<") < subFieldSizeY ("<< subFieldSizeX <<")"  << endl;
        cerr << "You can also use Aperture Photometry (AP) instead of Weighted Mask Photometry (WM) ."  << endl;
        exit(1);

    }
    
    //If 2*subPixelsPerPixel < subfield size
    if (photometryMethod == "WM" && ((subPixelsPerPixel+1)*2 > subFieldSizeX || (subPixelsPerPixel+1)*2 > subFieldSizeY) )
    {
        cerr << "\nError (PreProcessingPhotometry::PreProcessingPhotometry()): You must provide a higher subfield size or a smaller subPixel value because  "<< endl;
        cerr << "2*subPixelsPerPixel (2*"<< subPixelsPerPixel <<")  > subFieldSizeX ("<< subFieldSizeX <<") or "  
                "2*subPixelsPerPixel (2*"<< subPixelsPerPixel <<") > subFieldSizeY ("<< subFieldSizeX <<")"<< endl;
        cerr << "You can also use Aperture Photometry (AP) instead of Weighted Mask Photometry (WM) in the input parameters photometry file."  << endl;
        exit(1);

    }
    
    //Make directory to write the data. 
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();  
    photometryDirName = outputPath + "/" + prefix + "/Photometry" + prefix;
    photometryPlotsDir = photometryDirName + "/" + "PLOTS";
    
    LogManager::log<< "    Creating output directories";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //Create Photometry output directory.
    if (!FileUtilities::doMkdir(photometryDirName.c_str(), S_IRWXU))
    {
        cerr << "\nError (PreProcessingPhotometry::PreProcessingPhotometry()): Cannot create directory " << photometryDirName << endl;
        exit(1);
    }
    
    //Create  output directory for photometry plots.
	if (!FileUtilities::doMkdir(photometryPlotsDir.c_str(), S_IRWXU))
	{
        cerr << "\nError (PreProcessingPhotometry::PreProcessingPhotometry()): Cannot create directory " << photometryPlotsDir << endl;
		exit(1);
	}

    //Setting photometry parameters
	LogManager::log << "    Setting photometry parameters" ;
    GlobalVariables::logManager.LogManagerShowLog();    
    
    //Remove all photometry files in directory
	vector<string> starFiles;
	string dummy;
	FileUtilities::getDir(photometryDirName, starFiles, ".phot");
    
	for (uint i = 0; i < starFiles.size(); i++)
	{
		dummy = photometryDirName + "/" + starFiles[i];
		remove(dummy.c_str());
	}
    
	dummy = photometryDirName + "/star_coordinates.dat";
	remove(dummy.c_str());
    
    
    LogManager::log<< "    Output directory of photometry "<< photometryDirName << endl;
    LogManager::log<< "    Output directory of plots  "<< photometryPlotsDir;
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

        
    //Set PreProcessing photometry parameters into the PhotometryDataSet
    datasetPhotometry.datasetSetPhotometryParams(photometryDirName, photometryPlotsDir);
    
}
//==============================================================================








