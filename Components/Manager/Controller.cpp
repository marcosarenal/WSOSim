///////////////////////////////////////////////////////////
//  Controller.cpp
//  Implementation of the Class Controller
//  Created on:      23-Oct-2012 1:59:58 PM
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



#include "Controller.h"



//==============================================================================
/**
 * Constructor method
 */
Controller::Controller(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
Controller::~Controller(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================

/**
 * This is the main function of the Controller class. Depending on the user
 * input arguments the runController() chooses which kind of processing must
 * be launched by calling different functions.
 *
 * The options are:
 * 1) CCD Simulation
 * 2) CMOS Simulation
 * 3) Photometry Simulation
 * 4) CCD and Photometry Simulation
 * 5) CMOS and Photometry Simulation
 * 6) WUVS simulations
 *
 */
void Controller::runController(int argc, char ** argv)
{
    //Showing Presentation
    GlobalVariables::logManager.LogManagerPresentation();
    
    string parameterFile = string(argv[2]);
    
    
    //This option calls the CCD and Photometry processing.
    if (argc == 5 && string(argv[1]) == "-s" && string(argv[3]) == "-p")
    {
    //Identifying photometry parameter file.
    string photometryParameterFile = string(argv[4]);

    LogManager::log << "Starting CCD Simulation and Photometry";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runCCDandPhotometryController(parameterFile, photometryParameterFile);
    }
    

    
    //This option calls the CMOS and Photometry processing.
    else if (argc == 5 && string(argv[1]) == "-m" && string(argv[3]) == "-p")
    {
    //Identifying photometry parameter file.
    string photometryParameterFile = string(argv[4]);

    LogManager::log << "Starting CMOS Simulation and Photometry";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runCMOSandPhotometryController(parameterFile, photometryParameterFile);
    }
    
    //This option calls the CCD processing.
    else if (argc == 3 && string(argv[1]) == "-s")
    {
    LogManager::log << "Starting CCD Simulation";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runCCDController(parameterFile);
    }
    
    
    //This option calls the CMOS processing.
    else if (argc == 3 && string(argv[1]) == "-m")
    {
    LogManager::log << "Starting CMOS Simulation";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runCMOSController(parameterFile);
    }
    
    //This option calls the Photometry processing.
    else if (argc == 4 && string(argv[1]) == "-p")
    {
    //Identifying photometry parameter file.
    string photometryParameterFile = string(argv[2]);

    LogManager::log << "Starting Photometry Simulation";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runPhotometryController(parameterFile, photometryParameterFile);
    }

 
    //This option calls the WUVS processing.
    else if (argc == 3 && string(argv[1]) == "-w")
    {
    LogManager::log << "Starting WUVS Simulation";
    GlobalVariables::logManager.LogManagerShowLog();
    Controller::runWUVSController(parameterFile);
    }
 

    
    //In case of incorrect input parameters usage explanation is shown to user.
    else
    {
            cerr << "Usage: " << argv[0] << " [options] [parameter file] [[ccd input file]]" << endl;
            cerr << "options: " << endl;
            cerr << "-s ... Run CCD simulation" << endl;
            //cerr << "-m ... Run CMOS simulation" << endl;
            cerr << "-p ... Compute photometry" << endl;
            cerr << "parameter file: An xml-file with input parameters" << endl;
            cerr << "photometry parameter file: An XML-file with photometry input parameters."<< endl;
            cerr << endl << "Examples:" << endl;
            cerr << "Run a CCD simulation: ./PLATOsim -s /home/PLATOSim/inputfiles/ccd_parameters.xml" << endl;
            cerr << "Run a CCD simulation and make photometry: ./PLATOSim -s /home/PLATOSim/inputfiles/ccd_parameters.xml -p /home/PLATOSim/inputfiles/photometry_parameters.xml" << endl;
    }
    
    return;
}
//==============================================================================




//==============================================================================
/**
 * This function is called to perform the complete CCD and Photometry processing
 * and is in charge of trigger the PreProcessing, Processing and PostProcessing.
 */
void Controller::runCCDandPhotometryController(string parameterFile, string photometryParameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    //PHOTOMETRY DATASET DECLARATION
    DataSetPhotometry datasetPhotometry;
        
    //PRE-PROCESSING
    LogManager::log << "Starting PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //PreProcessing Common initialization for reading input parameter file and create output dir.
    PreProcessingCommon m_PreProcessingCommon;
    
    m_PreProcessingCommon.PreProcessingCommonCCD(m_DataSet, parameterFile);
    m_PreProcessingCommon.PreProcessingCommonPhotometry(datasetPhotometry, photometryParameterFile);
	
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
        
    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();

    PreProcessingCCD m_PreProcessingCCD(m_DataSet);
    PreProcessingPSF m_PreProcessingPSF(m_DataSet);
    LogManager::log << "End PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
        
    //PREPROCESSING PHOTOMETRY
    LogManager::log << "Starting Photometry PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    //Note: PreProcessing Photometry is performed BEFORE Processing for not to wait until the end of the Processing for checking wrong parameters.        
    PreProcessingPhotometry preProcessingPhotometry(m_DataSet, datasetPhotometry, photometryParameterFile);
    
    ProcessingPhotometry processingPhotometry; // Commented when processing photometry on-the-go

    
    //PROCESSING
    LogManager::log << "Starting CCD Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    ProcessingCCD m_ProcessingCCD;          //Declaration of the ProcessingCCD class.
    m_ProcessingCCD.processingCCDPipeline(m_DataSet); //Perform photometry from FITS files
//    m_ProcessingCCD.processingCCDPipelineOnTheGoPhotom(m_DataSet, datasetPhotometry); //Perform photometry on-the-go
    
    LogManager::log << "End CCD Processing";
    GlobalVariables::logManager.LogManagerShowLog();
           

    //PROCESSING PHOTOMETRY
    LogManager::log << "Starting Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    processingPhotometry.processingPhotometryPipeline(m_DataSet, datasetPhotometry); //Perform photometry from FITS files
    
    LogManager::log << "End Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //POST-PROCESSING
    LogManager::log << "Starting PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    PostProcessing m_PostProcessing(m_DataSet, datasetPhotometry);
    m_PostProcessing.postProcessingGenerateNoisePlots();
    
    LogManager::log << "End PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================




//==============================================================================
/**
 * This function is called to perform the complete CMOS and Photometry processing
 * and is in charge of trigger the PreProcessing, Processing and PostProcessing.
 */
void Controller::runCMOSandPhotometryController(string parameterFile, string photometryParameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    //PHOTOMETRY DATASET DECLARATION
    DataSetPhotometry datasetPhotometry;
    
    //PRE-PROCESSING
    LogManager::log << "Starting PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //PreProcessing Common initialization for reading input parameter file and create output dir.
    PreProcessingCommon m_PreProcessingCommon;
    
    //        m_PreProcessingCommon.PreProcessingCommonCMOS(m_DataSet, parameterFile);
    m_PreProcessingCommon.PreProcessingCommonPhotometry(datasetPhotometry, photometryParameterFile);
    
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
    
    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();

    PreProcessingCMOS m_PreProcessingCMOS;
    PreProcessingPSF m_PreProcessingPSF(m_DataSet);
    LogManager::log << "End PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //PREPROCESSING PHOTOMETRY
    LogManager::log << "Starting Photometry PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PreProcessingPhotometry preProcessingPhotometry(m_DataSet, datasetPhotometry, photometryParameterFile);
    //Note: PreProcessing Photometry is performed BEFORE Processing for not to wait until the end of the Processing for checking wrong parameters.
    
    
    //PROCESSING
    LogManager::log << "Starting CMOS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    ProcessingCMOS m_ProcessingCMOS;
    LogManager::log << "End CMOS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
       
    
    //POST-PROCESSING
    LogManager::log << "Starting PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PostProcessing m_PostProcessing(m_DataSet, datasetPhotometry);
    LogManager::log << "End PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();



    //PROCESSING PHOTOMETRY
    LogManager::log << "Starting Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    ProcessingPhotometry processingPhotometry;
    processingPhotometry.processingPhotometryPipeline(m_DataSet, datasetPhotometry);
    
    LogManager::log << "End Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();

    
}
//==============================================================================




//==============================================================================
/**
 * This function is called to perform the complete CCD processing and is in charge of
 * trigger the PreProcessing, Processing and PostProcessing.
 */
void Controller::runCCDController(string parameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    
    //PRE-PROCESSING
    LogManager::log << "Starting PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //PreProcessing Common initialization for reading input parameter file and create output dir.
    PreProcessingCommon m_PreProcessingCommon;
    m_PreProcessingCommon.PreProcessingCommonCCD(m_DataSet, parameterFile);
    
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
    
    PreProcessingCCD m_PreProcessingCCD(m_DataSet);
    PreProcessingPSF m_PreProcessingPSF(m_DataSet);
    LogManager::log << "End PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();

    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();
    
    
    //PROCESSING
    LogManager::log << "Starting CCD Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    ProcessingCCD m_ProcessingCCD;          //Declaration of the ProcessingCCD class.
    m_ProcessingCCD.processingCCDPipeline(m_DataSet);
    
    LogManager::log << "End CCD Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //POST-PROCESSING
    LogManager::log << "Starting PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PostProcessing m_PostProcessing();
    LogManager::log << "End PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================




//==============================================================================
/**
 * This function is called to perform the complete CMOS processing and is in charge of
 * trigger the PreProcessing, Processing and PostProcessing.
 */
void Controller::runCMOSController(string parameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    
    //PRE-PROCESSING
    LogManager::log << "Starting PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //PreProcessing Common initialization for reading input parameter file and create output dir.
    PreProcessingCommon m_PreProcessingCommon;
    //        m_PreProcessingCommon.PreProcessingCommonCMOS(m_DataSet, parameterFile);
    
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
    
    PreProcessingCMOS m_PreProcessingCMOS;
    PreProcessingPSF m_PreProcessingPSF(m_DataSet);
    LogManager::log << "End PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();
    
    //PROCESSING
    LogManager::log << "Starting CMOS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    ProcessingCMOS m_ProcessingCMOS;
    LogManager::log << "End CMOS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //POST-PROCESSING
    LogManager::log << "Starting PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PostProcessing m_PostProcessing();
    LogManager::log << "End PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================





//==============================================================================
/**
 * This function is called to perform the photometry processing and is in charge of
 * trigger the PreProcessing, Processing and PostProcessing.
 */
void Controller::runPhotometryController(string parameterFile, string photometryParameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    //PHOTOMETRY DATASET DECLARATION
    DataSetPhotometry datasetPhotometry;
    
    //PRE-PROCESSING
    LogManager::log << "Starting PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //PreProcessing Common initialization for reading input parameter file and create output dir.
    PreProcessingCommon m_PreProcessingCommon;
    m_PreProcessingCommon.PreProcessingCommonPhotometry(datasetPhotometry, photometryParameterFile);
    
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
    
    PreProcessingPSF m_PreProcessingPSF(m_DataSet);
    PreProcessingPhotometry preProcessingPhotometry(m_DataSet, datasetPhotometry, photometryParameterFile);
    LogManager::log << "End PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();
    
    //PROCESSING
    LogManager::log << "Starting Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    ProcessingPhotometry processingPhotometry;
    processingPhotometry.processingPhotometryPipeline(m_DataSet, datasetPhotometry);
    LogManager::log << "End Photometry Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //POST-PROCESSING
    LogManager::log << "Starting PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PostProcessing m_PostProcessing(m_DataSet, datasetPhotometry);
    LogManager::log << "End PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================





//==============================================================================
/**
 * This function is called to perform the simulation of the WUVS Spectrographs 
 * on board the WSO-UV spacecraft. This function takes as an input a sinthetic spectrograph 
 * image and is in charge of trigger its PreProcessing, Processing and PostProcessing.
 */
void Controller::runWUVSController(string parameterFile)
{
    //DATASET DECLARATION
    DataSet m_DataSet;
    
    //PRE-PROCESSING
    LogManager::log << "Starting WUVS PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    
     //PreProcessing Common initialization for reading the CCD input parameter file and create output dir.
    //TODO: The call to the PreProcessingCommonCCD for reading the CCD input parameter file must be replaced 
    //      for a new method to read a new WUVS input parameter file. Meanwhile, we will be using the old
    //      CCD parameter file.
    PreProcessingCommon m_PreProcessingCommon;
    m_PreProcessingCommon.PreProcessingCommonCCD(m_DataSet, parameterFile);
    
    //Creating the output file with the output path and name indicated into the xml input file.
    GlobalVariables::logManager.LogManagerGenerateLogFile(m_DataSet.datasetGetOutputPath() + "/" + m_DataSet.datasetGetPrefix() +"/"+ "simulation.info");
    
    PreProcessingWUVS m_PreProcessingWUVS(m_DataSet);
    LogManager::log << "End WUVS PreProcessing";
    GlobalVariables::logManager.LogManagerShowLog();

    //Showing Presentation    
    GlobalVariables::logManager.LogManagerPresentation();
    GlobalVariables::logManager.LogManagerAppendLog();
    
    
    //PROCESSING
    LogManager::log << "Starting WUVS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    ProcessingWUVS m_ProcessingWUVS;          //Declaration of the ProcessingWUVS class.
    m_ProcessingWUVS.processingWUVSPipeline(m_DataSet);
    
    LogManager::log << "End WUVS Processing";
    GlobalVariables::logManager.LogManagerShowLog();
    
    
    //POST-PROCESSING
    LogManager::log << "Starting WUVS PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
    PostProcessing m_PostProcessing();
    LogManager::log << "End WUVS PostProcessing";
    GlobalVariables::logManager.LogManagerShowLog();
 
}
//==============================================================================


