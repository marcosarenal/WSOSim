///////////////////////////////////////////////////////////
//  StepSaturation.cpp
//  Implementation of the Class StepSaturation
//  Created on:      23-Oct-2012 2:00:02 PM
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



#include "StepSaturation.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepSaturation::StepSaturation(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepSaturation::~StepSaturation(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the Full well saturation to pixelMap.
 */
void StepSaturation::StepSaturationfullWellapplication (DataSet &m_DataSet)
{
    
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	

    //Retrieving parameters from DataSet        
    fullWellSat = p_DataSet->datasetGetfullWellSat();

    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();

                  
   double e, diff;
   double sat = double(fullWellSat);
   int jmod;

   //compute transfer of electrons
   for ( int i = 0;i < pixelMap.rows();i++ )
   {
      for ( int j = 0;j < pixelMap.cols();j++ )
      {
         e = pixelMap ( i, j );
         if ( e > sat ) //well is saturated -> distribute the electrons in the wells above and below evenly until no saturation
         {
            //transfer down
            jmod = j;
            diff = ( e - sat ) / 2.; //the half down

            while ( diff > 0 && jmod < pixelMap.cols() )
            {
               pixelMap ( i, jmod ) -= diff;
               jmod++;
               if ( jmod < pixelMap.cols() )
               {
                  pixelMap ( i, jmod ) += diff;
                  if ( pixelMap ( i, jmod ) - diff < sat )
                  { 
                      diff = pixelMap ( i, jmod ) - sat;
                  }
               }
            }


            //transfer up
            jmod = j;
            e = pixelMap ( i, j );
            diff = e - sat; //the rest up
            while ( diff > 0 && jmod >= 0 )
            {
               pixelMap ( i, jmod ) -= diff;
               jmod--;
               if ( jmod >= 0 )
               {
                  pixelMap ( i, jmod ) += diff;
                  if ( pixelMap ( i, jmod ) - diff < sat ) 
                  {
                      diff = pixelMap ( i, jmod ) - sat;
                  }
               }
            }

         }
      }
   }
                        
       
   //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);
    
    
   LogManager::log <<"    Successfully added Full Well Saturation effect.";
   GlobalVariables::logManager.LogManagerShowLog();     

}




//==============================================================================
/**
 * This function applies the Digital saturation to pixelMap.
 */
void StepSaturation::StepSaturationdigitalapplication (DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;       	

    //Retrieving parameters from DataSet
    digitalSat = p_DataSet->datasetGetdigitalSat();     

    digitalSat = double(digitalSat);


    //Retrieving the BIAS register map from DataSet        
    //Initialize smearingMap
    biasRegisterMap.resize(m_DataSet.datasetGetbiasRegisterMap().extent(0), m_DataSet.datasetGetbiasRegisterMap().extent(1));
    biasRegisterMap = 0.0;
    biasRegisterMap = m_DataSet.datasetGetbiasRegisterMap();

    //Retrieving the smearing map from DataSet
    smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
    smearingMap = 0.0;
    smearingMap = m_DataSet.datasetGetsmearingMap();


    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();


    //Used Blitz function "where" 
    //Wherever pixelMap > digitalSat is true, digitalSat is returned;
    //where pixelMap > digitalSat is false, pixelMap is returned.
    pixelMap = where ( pixelMap > digitalSat, digitalSat, pixelMap);

    //Wherever biasRegisterMap > digitalSat is true, digitalSat is returned;
    //where biasRegisterMap > digitalSat is false, biasRegisterMap is returned.
    biasRegisterMap = where (biasRegisterMap > digitalSat, digitalSat, biasRegisterMap );

    //Wherever smearingMap > digitalSat is true, digitalSat is returned;
    //where smearingMap > digitalSat is false, smearingMap is returned.
    smearingMap = where ( smearingMap > digitalSat, digitalSat, smearingMap );


    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);

    //Setting the new calculated biasRegisterMap into the DataSet
    m_DataSet.datasetSetbiasRegisterMap(biasRegisterMap);

    //Setting the new calculated smearingMap into the DataSet
    m_DataSet.datasetSetsmearingMap(smearingMap);


    LogManager::log <<"    Successfully added Digital Saturation effect.";
    GlobalVariables::logManager.LogManagerShowLog();     


}
//==============================================================================