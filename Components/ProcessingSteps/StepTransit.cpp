///////////////////////////////////////////////////////////
//  StepTransit.cpp
//  Implementation of the Class StepTransit
//  Created on:      23-Sep-2016 11:38:59 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
 * This file is part of the WSO Simulator (WSOSim).
 *
 * WSOSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * WSOSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with WSOSim.  If not, see <http://www.gnu.org/licenses/>.
 */




#include "StepTransit.h"


//==============================================================================
/**
 * Constructor method
 */
StepTransit::StepTransit(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepTransit::~StepTransit(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the modification in the flux of the source when the transit of the exoplanet takes place.
 * The reduction in flux is calculated as the fraction of the surface of the star ocultated by the surface of the planet.
 * The transit times are retrieved from the array inTransitArray with the beguining and end of each transit as in the example:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...) 
 */
void StepTransit::StepTransitapplication(DataSet &m_DataSet, double startTime)
{   
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    hostStarTransitRA = p_DataSet->datasetGethostStarTransitRA();       
    hostStarTransitDec = p_DataSet->datasetGethostStarTransitDec();       
    hostStarRadius = p_DataSet->datasetGethostStarRadius();       
    exoplanetRadius = p_DataSet->datasetGetexoplanetRadius();       
    hostStarMagnitude = p_DataSet->datasetGethostStarMagnitude();       
    hostStarID = p_DataSet->datasetGethostStarID();       
    exoplanetOrbitalPeriod = p_DataSet->datasetGetexoplanetOrbitalPeriod();
      
    
    //Retrieving the inTransitArray from the DataSet
    inTransitArray.resize(p_DataSet->datasetGetinTransitArray().extent(0),p_DataSet->datasetGetinTransitArray().extent(1));
    inTransitArray = 0.0;
    inTransitArray = p_DataSet->datasetGetinTransitArray(); 

    inTransitArray = m_DataSet.datasetGetinTransitArray();
                 
             
    //Check whether the exposition startTime corresponds with a transit period on the inTransitArray
    
    //Check the aproximate startTime position in inTransitArray to iterate there instead of on the whole inTransitArray to reduce computing time.
    //Given the startTime and the exoplanet Orbital Period, we can aproximate the position of the startTime to verify whether this belongs to a transit or a non-transit period.
    orbitalPeriodInSeconds = exoplanetOrbitalPeriod * 24 * 60 * 60;
            
    // Retrieving the starCatalogue: contains for each star, its X and Y position in CCD, magnitude, RA,
    // declination and identification number.
    starCatalogue.resize(p_DataSet->datasetGetStarCatalogue().extent(0),p_DataSet->datasetGetStarCatalogue().extent(1));
    //starCatalogue = 0.0;
    starCatalogue = p_DataSet->datasetGetStarCatalogue();

     //Due to its structure inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...), init transit are in even positions of the array (that's why it's multiplied by 2)
    int firstIteration = floor(startTime/orbitalPeriodInSeconds) - 1;
 
    //Check that we are not in the first transit (startTime < orbitalPeriodInSeconds) 
    if(firstIteration<0)
    {
        firstIteration = 0;
    }
                
    //For the closer positions of startTime in the inTransitArray:
    for (int iterStartTime = firstIteration; iterStartTime < firstIteration + 2; iterStartTime++)
    {
        //If the startTime is in between an even and an odd position (and NOT between an odd and an even position), meaning that it is on transit
        if(inTransitArray(2*iterStartTime) < startTime && startTime < inTransitArray(2*iterStartTime+1))
        {
            //Calculate the aparent magnitude of the host star during the transit
            //F(transit)/F(no transit)=1-(exoplanetRadius^2/hostStarRadius^2)
            transitMagnitude = hostStarMagnitude - (1/0.4)*log10(1-(exoplanetRadius*exoplanetRadius/(hostStarRadius*hostStarRadius))); 


            //Modify the magnitude of the transit host star in the star catalogue so it takes into account the dimmering due to the exoplanet transit
            starCatalogue(hostStarID, 3) = transitMagnitude;
    
            LogManager::log <<"    This exposure will include an exoplanet transit:"<< std::endl;
            LogManager::log <<"    Host star ID in the Star Catalogue: "<< hostStarID +1<< std::endl;//+1 because the starCatalogue array starts in position 0.
            LogManager::log <<"    Magnitude during transit: "<< transitMagnitude<< std::endl;//+1 because the starCatalogue array starts in position 0.
            LogManager::log <<"    Transit start and end [s]: "<< inTransitArray(2*iterStartTime)<<"--> "<< inTransitArray(2*iterStartTime+1)<<std::endl;//+1 because the starCatalogue array starts in position 0.
            GlobalVariables::logManager.LogManagerAppendLog();

            //Set the starCatalogue into the DataSet            
            m_DataSet.datasetSetStarCatalogue(starCatalogue); 
        }
        else
        {   
            //Set the original hostStarMagnitude magnitude of the transit host star in the star catalogue 
            starCatalogue(hostStarID, 3) = hostStarMagnitude;

            //Set the starCatalogue into the DataSet            
            m_DataSet.datasetSetStarCatalogue(starCatalogue); 

            
        }
    }
        
   
    LogManager::log <<"    Successfully implemented exoplanet transit.";
    GlobalVariables::logManager.LogManagerShowLog();
}   