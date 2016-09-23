///////////////////////////////////////////////////////////
//  ParamsTransit.cpp
//  Implementation of the Class ParamsTransit
//  Created on:      21-Sep-2016 11:38:59 PM
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




#include "ParamsTransit.h"


//==============================================================================
/**
 * Constructor method
 */
ParamsTransit::ParamsTransit(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsTransit::~ParamsTransit(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * Sets the time when the transit of the exoplanet takes place.
 * It is set as the array inTransitArray with the beguining and end of each transit as in the example:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...) 
 */
void ParamsTransit::ParamsTransitsetting(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    //Retrieve transit parameters from dataset
    hostStarTransitRA = p_DataSet->datasetGethostStarTransitRA();
    hostStarTransitDec = p_DataSet->datasetGethostStarTransitDec();
    hostStarRadius = p_DataSet->datasetGethostStarRadius();
    exoplanetRadius = p_DataSet->datasetGetexoplanetRadius();
    exoplanetOrbitalPeriod = p_DataSet->datasetGetexoplanetOrbitalPeriod();
    planetaryOrbitSemiaxis = p_DataSet->datasetGetplanetaryOrbitSemiaxis();
    planetaryOrbitInclination = p_DataSet->datasetGetplanetaryOrbitInclination();
        
    //Check whether the input parameters are valid
    ParamsTransit::ParamsTransitcheckInputParameters(m_DataSet);

    //Calculate the beginning and end of each transit during the whole simulation time. 
    ParamsTransit::ParamsTransitcalculation(m_DataSet);
    
    //Setting the inTransitArray  (in s) calculated in the DataSet 
    p_DataSet->datasetSetinTransitArray(inTransitArray); 
//    //Set the calculated inTransitArray in the DataSet
//    m_DataSet.datasetSetinTransitArray(inTransitArray);	

}
//==============================================================================


//==============================================================================
/**
 * This function checks that every transit parameter is valid for performing the simulation.
 */
void ParamsTransit::ParamsTransitcheckInputParameters(DataSet &m_DataSet)
{
    performExoTransit = p_DataSet->datasetGetperformExoTransit();
	if (performExoTransit)
	{

            if (hostStarTransitRA < 0 || hostStarTransitRA > 360)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Right Ascension of the transit host star must be between 0 and 360." << endl;
		exit(1);
            }

            if (hostStarTransitDec < -90 || hostStarTransitDec > 90)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Declination of the transit host star must be between -90 and +90." << endl;
		exit(1);
            }
            
            if (hostStarRadius < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the host star must be > 0." << endl;
		exit(1);
            }
            
            if (exoplanetRadius < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the transit planet  must be > 0." << endl;
		exit(1);
            }
                        
            if (hostStarRadius < exoplanetRadius )
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Please, check that the radius of the transit planet is smaller than the radius of its host star." << endl;
		exit(1);
            }

            if (exoplanetOrbitalPeriod < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Exoplanet Orbital Period must be > 0." << endl;
		exit(1);
            }
            
            if (planetaryOrbitSemiaxis < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Semiaxis must be > 0." << endl;
		exit(1);
            }

            if (planetaryOrbitInclination < 0 || planetaryOrbitInclination > 180 )
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Inclination must be between 0 and 180." << endl;
		exit(1);
            }
            
            //Calculation of the impact parameter to check whether there will be a transit or not
            if ( (planetaryOrbitSemiaxis / ( 0.00465047 * hostStarRadius ) ) * std::cos(planetaryOrbitInclination * Constants::Pi / 180) > 1 )
            {
                cout<<"_________________________" << (planetaryOrbitSemiaxis / ( 0.00465047 * hostStarRadius ) ) * std::cos(planetaryOrbitInclination * Constants::Pi / 180)<<endl;
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Under these conditions there will not be an observable transit. Please modify them (an orbit inclination closer to 90 deg would help)." << endl;
		exit(1);
            }


	}
    
    
}
//==============================================================================

//==============================================================================
/**
 * Compute the beginning and end of each transit during the whole simulation time.
 * This method takes the exoplanet orbital parameters and sets in the array inTransitArray 
 * the beginning and end of each transit in the form:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...).
 * The beguinning of the transit starts in a random time.
 */
void ParamsTransit::ParamsTransitcalculation(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieve parameters from the DataSet
    exposureTime = m_DataSet.datasetGetexposureTime();
    readOutTime = m_DataSet.datasetGetreadOutTime();
    integrationTime = exposureTime + readOutTime;    
    numExposures = m_DataSet.datasetGetnumExposures();
       
    
    //Calculate the duration of the whole simulation time (in seconds)
    totalSimulationDuration = numExposures * integrationTime;

    //Transform the orbital period from days to seconds
    orbitalPeriodInSeconds = exoplanetOrbitalPeriod * 24 * 60 * 60;
    
    //Calculate the number of trasits in that amount of time
    numberOfTransits = (int)floor(totalSimulationDuration / orbitalPeriodInSeconds);
    
    cout<<"_________ ParamsTransit::ParamsTransitcalculation __________"<< numberOfTransits <<endl;
    
    //Check whether there will be at least one transit during the simulation duration    
    
    //Resize the transit array to twicw the number of transits in the whole simulation
    //in order to allocate beginning and end time of each of them
    inTransitArray.resize(2 * numberOfTransits);
    
    //If the numberOfTransits is zero, there might be a transit or not. The inTransitArray must be resized to 2 anyway. 
    if (numberOfTransits == 0)
    {
        inTransitArray.resize(2);        
    }
    

    //Seed the random generator
    std::srand(DataSet::seedRNG++);
        
    //Set the beginning of the 1st transit in a random time < orbital period
    inTransitArray(0) = ( std::rand() / (double)RAND_MAX) * orbitalPeriodInSeconds;

    //Calculate the transit duration as in https://www.paulanthonywilson.com/exoplanets/exoplanet-detection-techniques/the-exoplanet-transit-method/
    
    //Calculate planetaryOrbitSemiaxis in Solar radius units (1 solar radius = 0.00465047 AU)
    planetaryOrbitSemiaxisInSolarRadius = planetaryOrbitSemiaxis / 0.00465047;
    //Calculate the impact parameter
    impactParameter = (planetaryOrbitSemiaxisInSolarRadius /  hostStarRadius  ) * std::cos(planetaryOrbitInclination * Constants::Pi / 180);
    
    //Calculate the transit duration
    transitDuration = (orbitalPeriodInSeconds / Constants::Pi) * std::asin(std::sqrt(std::pow(hostStarRadius+exoplanetRadius,2)-std::pow(impactParameter * hostStarRadius,2))/planetaryOrbitSemiaxisInSolarRadius);
    
    //Set the time of the transit end on the inTransitArray as the beguinning of the transit plus the transit duration
    inTransitArray(1) = inTransitArray(0) + transitDuration;
    
    //The first transit is already set. 
    if (numberOfTransits > 1)
    {
        //For the following transits:
        for (int iterTransit = 1; iterTransit < numberOfTransits; iterTransit++)
        {
            cout<<"______ iterTransit _____"<< iterTransit <<endl;
            cout<<"______ inTransitArray(0) _____"<< inTransitArray(0) <<endl;
            cout<<"______ inTransitArray(1) _____"<< inTransitArray(1) <<endl;
            cout<<"______ orbitalPeriodInSeconds _____"<< orbitalPeriodInSeconds <<endl;
            inTransitArray(2*iterTransit) = inTransitArray(0) + iterTransit * orbitalPeriodInSeconds;
            inTransitArray(2*iterTransit + 1) = inTransitArray(1) + iterTransit * orbitalPeriodInSeconds;
                        
            
            cout<<"______ inTransitArray(2*iterTransit) _____"<< inTransitArray(2*iterTransit) <<endl;
            cout<<"______ inTransitArray(2*iterTransit+1) _____"<< inTransitArray(2*iterTransit+1) <<endl;

            
        }
    }
    
    cout<<"_________ ParamsTransit::ParamsTransitcalculation _____ transitDuration _____"<< transitDuration <<endl;
    cout<<"_________ ParamsTransit::ParamsTransitcalculation _____orbitalPeriodInSeconds_____"<< orbitalPeriodInSeconds <<endl;
    cout<<"_________ ParamsTransit::ParamsTransitcalculation _____inTransitArray(0)_____"<< inTransitArray(0) <<endl;
    cout<<"_________ ParamsTransit::ParamsTransitcalculation __inTransitArray________"<< inTransitArray <<endl;

   
  
}
//==============================================================================
