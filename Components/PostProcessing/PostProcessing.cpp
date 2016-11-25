///////////////////////////////////////////////////////////
//  PostProcessing.cpp
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



#include "PostProcessing.h"



//==============================================================================
/**
 * Constructor method
 */
PostProcessing::PostProcessing(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
PostProcessing::~PostProcessing(){}
//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
PostProcessing::PostProcessing(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    p_datasetPhotometry = &datasetPhotometry;


    //Initialize PSF map  
    psfMap.resize(m_DataSet.datasetGetPSFMap().extent(0),m_DataSet.datasetGetPSFMap().extent(1));
    psfMap = 0.0;
    //Retrieving the PSF Map from DataSet
    psfMap = m_DataSet.datasetGetPSFMap();

    //Retrieving parameters from DataSet
    exposureTime = p_DataSet->datasetGetexposureTime();
    flatfieldIntraPixelWidth = p_DataSet->datasetGetflatfieldIntraPixelWidth();
    //flatfieldIntraPixelWidth = 0;
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();
    outputDir = outputPath + "/" + prefix;    

    //Creating output FITS file
//        fileName = outputPath + "/" + prefix + "/" + prefix + mapName; 
//        
//        
//        
//    
//    //Initializing and retrieving the pixel map from DataSet
//    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
//    pixelMap = 0.0;
//    pixelMap = m_DataSet.datasetGetpixelMap();
//
//    //Initializing and retrieving the subpixel map from DataSet
//    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
//    subPixelMap = 0.0;
//    subPixelMap = m_DataSet.datasetGetsubPixelMap();
//    
//    //Initialize cteMap
//    cteMap.resize(m_DataSet.datasetGetCTEMap().extent(0), m_DataSet.datasetGetCTEMap().extent(1));
//    cteMap = 0.0;
//    //Retrieving the cteMap from DataSet
//    cteMap = m_DataSet.datasetGetCTEMap();

    //Initialize smearingMap
    smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
    smearingMap = 0.0;
    //Retrieving the smearingMap from DataSet
    smearingMap = m_DataSet.datasetGetsmearingMap();
//    
//        Trigger the  postProcessingWriteFITS
//        FileUtilities::FileUtilitiesWriteFITS(exposureTime, flatfieldIntraPixelWidth, pixelMap, outputDir + "/" + prefix + "PixelMap");
////      FileUtilities::FileUtilitiesWriteFITS(exposureTime, flatfieldIntraPixelWidth, cteMap, outputDir + "/" + prefix + "cteMap");
    FileUtilities::FileUtilitiesWriteFITS(exposureTime, flatfieldIntraPixelWidth, psfMap, outputDir + "/" + prefix + "psfMap");
////      FileUtilities::FileUtilitiesWriteFITS(exposureTime, flatfieldIntraPixelWidth, subPixelMap, outputDir + "/" + prefix + "subPixelMap");
    FileUtilities::FileUtilitiesWriteFITS(exposureTime, flatfieldIntraPixelWidth, smearingMap, outputDir + "/" + prefix + "smearingMap");


}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * 
 */
//void PostProcessing::postProcessingWriteCCDParams(){
//
//}
//==============================================================================




//==============================================================================
/**
 * 
 */
//void PostProcessing::postProcessingWriteCMOSParams(){
//
//}
//==============================================================================






//==============================================================================
/**
 * Wrapper class to handle XML files. Provides convenience methods to read and
 * write XML files.
 */
void PostProcessing::postProcessingWriteXML(){}
//==============================================================================





//==============================================================================
/**
 * Function to read the generated photometry files and generate noise statistics 
 * and plots. Each star photometry file contains in columns:
 * star*.phot ............. Lightcurve of each star
                                      Column description:
                                      1: time [s]
                                      2: Input magnitude
                                      3: Input flux
                                      4: Measured magnitude
                                      5: Measured normalized magnitude
                                      6: Measured flux
                                      7: Measured normalized flux
                                      8: Background flux
 * 
 */
void PostProcessing::postProcessingGenerateNoisePlots()
{
    
    photometryPlotsDir = p_datasetPhotometry->datasetGetphotometryPlotsDir();
    prefix = p_DataSet->datasetGetPrefix();	
    photometryDirName = p_datasetPhotometry->datasetGetphotometryDirName();
    exposureTime = p_DataSet->datasetGetexposureTime();
    background =  p_datasetPhotometry->datasetGetphotometryBackground();
    photometryNumTelescopes = p_datasetPhotometry->datasetGetphotometryNumTelescopes();

    std::string noise_statistics = photometryPlotsDir + "/noise_statistics.py";
	ofstream noiseScript(noise_statistics.c_str());
	if (!noiseScript.is_open())
	{
		std::cerr << "\nError (PostProcessing::postProcessingGenerateNoisePlots()): Unable to open output gnuscript file " << noise_statistics;
		exit(1);
	}
    
    LogManager::log << "    Creating script for python (" << noise_statistics << ")";
    GlobalVariables::logManager.LogManagerShowLog(); 
    
	noiseScript.precision(15);
	noiseScript << "#! /usr/bin/env python " << std::endl;
	noiseScript << "# -*- coding: utf-8 -*- " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "# This python scrip is automatically generated in the execution of the PLATOSim simulator.  " << std::endl;
	noiseScript << "#  " << std::endl;
	noiseScript << "# Every exposure generates a FITS file from which is measured the flux of every source;  " << std::endl;
	noiseScript << "# For every detected star in those FITS images, a photometry file ('.phot') is generated.  " << std::endl;
	noiseScript << "# Those photometry files contain the measured parameters in the FITS images, being added a new line  " << std::endl;
	noiseScript << "# for each of those images until the time-serie is completed with all the expositions. Therefore,  " << std::endl;
	noiseScript << "# each photometry file contains the measured parameters of each star for each exposure in the time-serie.  " << std::endl;
	noiseScript << "# Here are taken those input parameters in the photometry files (*.phot) from the "<< photometryDirName<<" folder " << std::endl;
	noiseScript << "# to generate two plots. The left hand-side plot shows the Measured magnitude of every star as a function of  " << std::endl;
	noiseScript << "# its Input magnitude. The Measured magnitude for each star is the mean value of all the measured magnitudes  " << std::endl;
	noiseScript << "# in the time-serie. " << std::endl;
	noiseScript << "# The right-hand-side plot shows the noise in ppm/sqrt(hour) as a function of the input magnitude. " << std::endl;
	noiseScript << "# The noise of each star is obtained using the standard deviation of the measured magnitude, normalized  " << std::endl;
	noiseScript << "# in the exposure time. " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "import pylab " << std::endl;
	noiseScript << "import numpy as np " << std::endl;
	noiseScript << "import matplotlib.pyplot as plt " << std::endl;
	noiseScript << "import glob " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "#name of the output folder containing the .phot files " << std::endl;
	noiseScript << "name='" << prefix << "' " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "#Definitions " << std::endl;
	noiseScript << "time=[]  " << std::endl;
	noiseScript << "input_mag=[]                 " << std::endl;
	noiseScript << "input_flux=[]                " << std::endl;
	noiseScript << "measured_mag=[]              " << std::endl;
	noiseScript << "#norm_measured_mag=[]         " << std::endl;
	noiseScript << "measured_flux=[]  " << std::endl;
	noiseScript << "norm_measured_flux=[]   " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "background_flux=" << background << std::endl;
//	noiseScript << "expTimeFactor = 1.0e6 / np.sqrt(3600. / " << exposureTime << ") # Factor to normalize measured noise ppm per hour " << std::endl;
	noiseScript << "expTimeFactor = 1.0e6 / np.sqrt(3600. / " << exposureTime << " * "<<photometryNumTelescopes<<") # Factor to normalize measured noise ppm per hour and number of telescopes" << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "#List of the photometry files (.phot) " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "file_list= glob.glob('" << photometryDirName << "/*.phot') " << std::endl;
	noiseScript << "mean_measured_mag=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "medianFluxInput=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "medianMagInput=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "medianFluxNorm=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "stdDevMag=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "medianFlux=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "stdDevFlux=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "#medianMag=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "#medianMagNorm=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "#stdDevInputMag=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "normFactor=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "#rmsMag=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "theoreticalPhotonNoise=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "ppmphObs=np.zeros(len(file_list)) " << std::endl;
	noiseScript << "ppmphPN=np.zeros(len(file_list)) " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "#Read each photometry file and set their parameters in arrays " << std::endl;
	noiseScript << "for i in range(len(file_list)): " << std::endl;
	noiseScript << "    f = open(file_list[i], 'r')   " << std::endl;
	noiseScript << "    header = f.readline()  " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << " " << std::endl;
	noiseScript << "#star*.phot ............. Lightcurve of each star " << std::endl;
	noiseScript << "#                                      Column description: " << std::endl;
	noiseScript << "#                                      1: time [s] " << std::endl;
	noiseScript << "#                                      2: Input magnitude " << std::endl;
	noiseScript << "#                                      3: Input flux " << std::endl;
	noiseScript << "#                                      4: Measured magnitude " << std::endl;
	noiseScript << "#                                      5: Measured normalized magnitude " << std::endl;
	noiseScript << "#                                      6: Measured flux " << std::endl;
	noiseScript << "#                                      7: Measured normalized flux " << std::endl;
	noiseScript << "#                                      8: Background flux " << std::endl;
	noiseScript << "              " << std::endl;
	noiseScript << "    #for each line in the .phot file " << std::endl;
	noiseScript << "    input_mag=[]              " << std::endl;
	noiseScript << "    input_flux=[]                " << std::endl;
	noiseScript << "    measured_mag=[]              " << std::endl;
	noiseScript << "    #norm_measured_mag=[]         " << std::endl;
	noiseScript << "    measured_flux=[] " << std::endl;
	noiseScript << "    norm_measured_flux=[]   " << std::endl;
	noiseScript << "     " << std::endl;
	noiseScript << "     " << std::endl;
	noiseScript << "     " << std::endl;	
    noiseScript << "    for line in f.xreadlines():" << std::endl;
	noiseScript << "        line = line.strip() " << std::endl;
	noiseScript << "        columns = line.split() " << std::endl;
	noiseScript << "         " << std::endl;
	noiseScript << "        time=(columns[0])      " << std::endl;
	noiseScript << "        input_mag.append(float(columns[1]))      " << std::endl;
	noiseScript << "        input_flux.append(float(columns[2]))      " << std::endl;
	noiseScript << "        measured_mag.append(float(columns[3]))      " << std::endl;
	noiseScript << "        #norm_measured_mag.append(float(columns[4]))      " << std::endl;
	noiseScript << "        measured_flux.append(float(columns[5]))     " << std::endl;
	noiseScript << "        norm_measured_flux.append(float(columns[6]))  " << std::endl;
	noiseScript << "         " << std::endl;
	noiseScript << "        mean_measured_mag[i]=np.mean(measured_mag) " << std::endl;
	noiseScript << "        medianMagInput[i]=np.mean(input_mag) " << std::endl;
	noiseScript << "        medianFluxInput[i]=np.mean(input_flux) " << std::endl;
	noiseScript << "        medianFluxNorm[i]=np.mean(norm_measured_flux) " << std::endl;
	noiseScript << "        #stdDevMag[i]=np.std(measured_mag) " << std::endl;
	noiseScript << "        medianFlux[i]=np.median(measured_flux) " << std::endl;
	noiseScript << "        stdDevFlux[i]=np.std(measured_flux) " << std::endl;
	noiseScript << "        #medianMag[i]=np.median(measured_mag) " << std::endl;
	noiseScript << "        #medianMagNorm[i]=np.median(norm_measured_mag) " << std::endl;
	noiseScript << "        #stdDevInputMag[i]=np.std(input_mag) " << std::endl;
	noiseScript << "        normFactor[i] = medianFluxNorm[i] / medianFlux[i] " << std::endl;
 	noiseScript << "	#rmsMag[i] = np.sqrt((medianMag[i] * medianMag[i]) + (stdDevMag[i]*stdDevMag[i])) " << std::endl;
	noiseScript << "	theoreticalPhotonNoise[i] = np.sqrt(medianFlux[i]) " << std::endl;
	noiseScript << "	ppmphPN[i] = (1. / theoreticalPhotonNoise[i]) * expTimeFactor " << std::endl;
    noiseScript << "	ppmphObs[i] = (stdDevFlux[i] * normFactor[i] / medianFluxInput[i]) * expTimeFactor #//measured noise ppm per hour (stdDevFlux/medianFlux) " << std::endl;
    noiseScript << "	 " << std::endl;
    noiseScript << "       " << std::endl;
    noiseScript << "    f.close() " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "  " << std::endl;
    noiseScript << "plt.figure(1, figsize=(20, 6)) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "ax1 = plt.subplot(1,2,1) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "ax1.plot(medianMagInput, mean_measured_mag, marker='+', color='red', markersize=4, linestyle='none') " << std::endl;
    noiseScript << "ax1.plot([8.5,14.7], [8.5,14.7], color='green', lw=0.5) " << std::endl;
    noiseScript << "#ax1.axis([6,15.1,6,15]) " << std::endl;
    noiseScript << "#ax1.set_xticks(np.arange(6,16)) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "plt.xlabel('Input Magnitude') " << std::endl;
    noiseScript << "plt.ylabel('Measured Magnitude') " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "## legend " << std::endl;
    noiseScript << "plt.legend((r'Observed magnitude', r'Input magnitude'), shadow = True , loc = 2) " << std::endl;
    noiseScript << "ltext = plt.gca().get_legend().get_texts() " << std::endl;
    noiseScript << "plt.setp(ltext[0], fontsize = 8, color = 'r') " << std::endl;
    noiseScript << "plt.setp(ltext[1], fontsize = 8, color = 'g') " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "plt.title(name,fontsize=14) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "ax2 = plt.subplot(1,2,2) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "ax2.plot(medianMagInput[:], ppmphObs[:],marker='+', color='red', markersize=4, linestyle='none') " << std::endl;
    noiseScript << "ax2.plot(medianMagInput[:], ppmphPN[:], marker=',',linestyle='none', color='green',lw=0.8) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "#axix " << std::endl;
    noiseScript << "#ax2.axis([8,15,10,10e4]) " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "# log scales: " << std::endl;
    noiseScript << "ax2.set_yscale('log') " << std::endl;
    noiseScript << "ax2.set_xscale('linear') " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "plt.xlabel('Input Magnitude') " << std::endl;
    noiseScript << "plt.ylabel('Measured noise (ppm/sqrt(hour))') " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "## legend " << std::endl;
    noiseScript << "plt.legend((r'Observed noise', r'theoretical photon noise'), shadow = True , loc = 2) " << std::endl;
    noiseScript << "ltext = plt.gca().get_legend().get_texts() " << std::endl;
    noiseScript << "plt.setp(ltext[0], fontsize = 8, color = 'r') " << std::endl;
    noiseScript << "plt.setp(ltext[1], fontsize = 8, color = 'g') " << std::endl;

    noiseScript << " " << std::endl;
    noiseScript << " " << std::endl;
    noiseScript << "plt.show() " << std::endl;
    noiseScript << "plt.savefig('" << photometryPlotsDir << "/noiseVSinput.pdf') " << std::endl;
    
    
     
	noiseScript.close();
    
	//execute the noise_statistics script
	std::string command = "chmod a+x " + noise_statistics + "; " + noise_statistics;
    
//	LogManager::log<< "    Noise plot generated in " << photometryPlotsDir ;
//    GlobalVariables::logManager.LogManagerAppendLogAndShow();  
    
	//system(command.c_str());
}
//==============================================================================

