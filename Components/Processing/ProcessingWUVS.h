///////////////////////////////////////////////////////////
//  ProcessingWUVS.h
//  Implementation of the Class ProcessingWUVS
//  Created on:      29-May-2015 13:00:00 PM
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



#ifndef PROCESSINGWUVS_H_
#define PROCESSINGWUVS_H_

#include "LogManager.h"
#include "DataSet.h"

#include "StepReadOutNoise.h"
#include "StepFlatField.h"
#include "StepRebin.h"
#include "StepBackground.h"
#include "StepCTE.h"
#include "StepChargeTransferSmearing.h"
#include "StepCosmics.h"
#include "StepElectronicOffset.h"
#include "StepGain.h"
#include "StepJitter.h"
#include "StepConvolvePSF.h"
#include "StepPhotonNoise.h"
#include "StepSaturation.h"
#include "StepPhotometryOnTheGo.h"
#include "PostProcessing.h"

/**
 * Class which calls the main methods for adding noise to the WUVS image.
 */
class ProcessingWUVS
{
public:
    
    ProcessingWUVS();
    virtual ~ProcessingWUVS();
    void processingWUVSPipeline(DataSet &m_DataSet);
    void processingWUVSPipelineOnTheGoPhotom(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry);
    
private:
    
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
    StepBackground m_StepBackground;               //Declaration of the StepBackground method.
    StepChargeTransferSmearing m_StepChargeTransferSmearing;        //Declaration of the StepChargeTransferSmearing method.
    StepCosmics m_StepCosmics;                      //Declaration of the StepCosmics method.
    StepCTE m_StepCTE;                              //Declaration of the StepCTE method.
    StepElectronicOffset m_StepElectronicOffset;    //Declaration of the StepElectronicOffset method.
    StepFlatField m_StepFlatField;                  //Declaration of the StepFlatField method.
    StepRebin m_StepRebin;                          //Declaration of the StepRebin method.
    StepGain m_StepGain;                            //Declaration of the StepGain method.
    StepJitter m_StepJitter;                        //Declaration of the StepJitter method.
    StepPhotonNoise m_StepPhotonNoise;              //Declaration of the StepPhotonNoise method.
    StepConvolvePSF m_StepConvolvePSF;              //Declaration of the StepConvolvePSF method.
    StepReadOutNoise m_StepReadOutNoise;            //Declaration of the StepReadOutNoise method.
    StepSaturation m_StepSaturation;                //Declaration of the StepSaturation method.    
    StepPhotometryOnTheGo m_StepPhotometryOnTheGo;  //Declaration of the StepSaturation method.

    
    bool   useJitter;                      //Parameter retrieved from DataSet.
    double flatfieldPixelNoise;            //Parameter retrieved from DataSet.
    double flatfieldWhiteNoise;            //Parameter retrieved from DataSet.
    double flatfieldIntraPixelWidth;       //Parameter retrieved from DataSet.
    int    numExposures;                   //Parameter retrieved from DataSet.
    string outputPath;                     //Parameter retrieved from DataSet.
    string prefix;                         //Parameter retrieved from DataSet.
    string outputDir;                      //Output directory path for the FITS files to be written.
    int    edgePixels;                     //Parameter retrieved from DataSet.
    int    numSmearingOverscanRows;        //Parameter retrieved from DataSet.
    int    numPrescanRows;                 //Parameter retrieved from DataSet.    

 
    double exposureTime;                   //Parameter retrieved from DataSet.
    double readOutTime;                    //Parameter retrieved from DataSet.
    double integrationTime;                // integrationTime = exposureTime + readOutTime.
    double startTime;                      //Starting time for each exposure to be applied the Jitter
     
    string log;                            //log string to be send to the LogManager.
     
    Array<float, 2>  pixelMap;                     //Blitz map set into the DataSet.        
    Array<float, 2>  subPixelMap;                  //Blitz map set into the DataSet.
    Array<float, 2>  initPixelMap;                 //Blitz map set into the DataSet.        
    Array<float, 2>  initSubPixelMap;              //Blitz map set into the DataSet.

    Array<string, 1> exposuresNamesArray;          //Array containing the exposures names for the FITS files.
    
    
};
#endif /* PROCESSINGWUVS_H_ */
