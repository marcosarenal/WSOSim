///////////////////////////////////////////////////////////
//  PreProcessingCommon.cpp
//  Implementation of the Class PreProcessingCommon
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




#include "PreProcessingCommon.h"



//==============================================================================
/**
 * Constructor method
 */
PreProcessingCommon::PreProcessingCommon(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
PreProcessingCommon::~PreProcessingCommon(){}
//==============================================================================




//==============================================================================
/**
 * The PreProcessingCommonCCD is in charge of reading the CCD input file and creating the output
 * directory using the path and name set in this input file; Those functions are common to 
 * any kind of processing and needs to be executed before any other.
 * Here is also copied the input XML file to the output directory.
 * @param m_DataSet DataSet object designed to store all the input parameters (among others)
 *        included into the xml parameters input file.
 * @param parameterFile Name of the xml parameters input file.
 */
void PreProcessingCommon::PreProcessingCommonCCD(DataSet &m_DataSet, string parameterFile)
{
    //Point to the DataSet
    p_DataSet = &m_DataSet;
            
    //Read parameters file.
    p_DataSet->datasetReadParameterFile(parameterFile);
    
    //Make directory to write the data. 
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();
    outputDir = outputPath + "/" + prefix;    
    
 
//    //Check if the output directory exists
//    if (FileUtilities::dirExists(outputDir.c_str()))    
//    {
//        string input = "";
//        LogManager::log << "\nWARNING: Output directory already exist; Continue will delete all its content! Continue? (y/n)" << endl;
//        GlobalVariables::logManager.LogManagerShowLog();
//        getline(cin, input);
//            
//        //If user says something different than yes, leave the simulation 
//        if (input != "y" && input != "Y")
//        {
//            cerr << "\nSimulation terminated "  << endl;
//            exit(1);
//        }
//        //if says yes, clean the output directory
//        else
//        {
//            string command = "rm -rf " +  outputDir;
//            int ret = system(command.c_str());
//            if (ret != 0)
//            {
//                cerr << "\nError (PreProcessingCommon::PreProcessingCommonCCD()): Cannot remove directory " << outputDir << endl;
//                exit(1);
//            }
//            
//        }
//    }

    
    //Creating output directory 
    if (!FileUtilities::doMkdir(outputDir.c_str(), S_IRWXU))
    {
        cerr << "\nError (PreProcessingCommon::PreProcessingCommonCCD()): Cannot create directory " << outputDir << endl;
        exit(1);
    }
    
    //Copy input XML file to output directory
    string command = "cp " + parameterFile + " " + outputDir + "/" + prefix + ".xml";
             
    int ret = system(command.c_str());
    if (ret != 0)
    {
        cerr << "\nError (PreProcessingCommon::PreProcessingCommonCCD()): Cannot copy directory " << outputDir << endl;
        exit(1);
    }            
            
}
//==============================================================================




//==============================================================================
/**
 * The PreProcessingCommonPhotometry is in charge of reading the input photometry file and creating the output
 * directory using the path and name set in this input file.
 * @param datasetPhotometry DataSetPhotometry object designed to store all the photometry 
 *      input parameters (among others) included into the xml photometry parameters input file.
 * @param photometryParameterFile Name of the xml photometry parameters input file.
 */
void PreProcessingCommon::PreProcessingCommonPhotometry(DataSetPhotometry &datasetPhotometry, string photometryParameterFile)
{

    //Read parameters file.
    datasetPhotometry.datasetPhotometryReadParameterFile(photometryParameterFile);
    

    
}
//==============================================================================