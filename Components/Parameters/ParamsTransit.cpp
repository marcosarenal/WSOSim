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
    ParamsTransit::ParamsTransittimeCalculation(m_DataSet);
        

    //Search for the transit host star in the star Catalogue
    ParamsTransit::ParamsTransitsearchHostStarOnCatalogue(m_DataSet);
    
    //Calculate the Flux of the transit host star 
    //ParamsTransit::ParamsTransitfluxCalculation(m_DataSet);
    
    //Setting the inTransitArray  (in s) calculated in the DataSet 
    p_DataSet->datasetSetinTransitArray(inTransitArray); 
        

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
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Right Ascension of the transit host star must be between 0 and 360." << std::endl;
		exit(1);
            }

            if (hostStarTransitDec < -90 || hostStarTransitDec > 90)
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Declination of the transit host star must be between -90 and +90." << std::endl;
		exit(1);
            }
            
            if (hostStarRadius < 0)
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the host star must be > 0." << std::endl;
		exit(1);
            }
            
            if (exoplanetRadius < 0)
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the transit planet  must be > 0." << std::endl;
		exit(1);
            }
                        
            if (hostStarRadius < exoplanetRadius )
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Please, check that the radius of the transit planet is smaller than the radius of its host star." << std::endl;
		exit(1);
            }

            if (exoplanetOrbitalPeriod < 0)
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Exoplanet Orbital Period must be > 0." << std::endl;
		exit(1);
            }
            
            if (planetaryOrbitSemiaxis < 0)
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Semiaxis must be > 0." << std::endl;
		exit(1);
            }

            if (planetaryOrbitInclination < 0 || planetaryOrbitInclination > 180 )
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Inclination must be between 0 and 180." << std::endl;
		exit(1);
            }
            
            //Calculation of the impact parameter to check whether there will be a transit or not
            if ( (planetaryOrbitSemiaxis / ( 0.00465047 * hostStarRadius ) ) * std::cos(planetaryOrbitInclination * Constants::Pi / 180) > 1 )
            {
		std::cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Under these conditions there will not be an observable transit. Please modify them (an orbit inclination closer to 90 deg would help)." << std::endl;
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
void ParamsTransit::ParamsTransittimeCalculation(DataSet &m_DataSet)
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
            //Set the beginning of each transit as the beginning of the 1st transit times the orbital period (even positions in the inTransitArray)
            inTransitArray(2*iterTransit) = inTransitArray(0) + iterTransit * orbitalPeriodInSeconds;
            //Set the end of each transit as the end of the 1st transit times the orbital period (odd positions in the inTransitArray)
            inTransitArray(2*iterTransit + 1) = inTransitArray(1) + iterTransit * orbitalPeriodInSeconds;
      
        }
    }
}

//==============================================================================

//==============================================================================
/**
 * This method searchs for the transit host star on the Star Catalogue based on 
 * the hostStarTransitRA and hostStarTransitDec parameters.
 * The method retrieves the position of the host star on the Star Catalogue and 
 * the magnitude that is provided there.
 * The hostStarTransitRA and hostStarTransitDec parameters must match with the 
 * Right Ascension and Declination on the Star Catalogue up 0.00001 degrees.
 */
void ParamsTransit::ParamsTransitsearchHostStarOnCatalogue(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieve parameters from the DataSet

    // Retrieving the starCatalogue: contains for each star, its X and Y position in CCD, magnitude, RA,
    // declination and identification number.
    starCatalogue.resize(p_DataSet->datasetGetStarCatalogue().extent(0),p_DataSet->datasetGetStarCatalogue().extent(1));
    starCatalogue = 0.0;
    starCatalogue = p_DataSet->datasetGetStarCatalogue();


    //Determine the number of sources on the Star Catalogue
    int NumStars = starCatalogue.rows();
        
    //Initialize parameter to check if the host star is found on the Star Catalogue
    bool transitHostStarMatch = false;
            
    //Initialize the star iterator in the Star Catalogue
    int starIterator = 0;
    
    //Iterate in the whole Star Catalogue until found a transit Host Star Match or end of the StarCatalogue
    while ( transitHostStarMatch == false && starIterator < NumStars)
    {
           
        //Check whether is a match of the hostStarTransitRA and hostStarTransitDec on the Star Catalogue and has not been found yet
        if (!transitHostStarMatch && fabs(starCatalogue(starIterator, 1) - hostStarTransitRA) < 0.00001 && fabs(starCatalogue(starIterator, 2) - hostStarTransitDec < 0.00001))
        {
            //Indicate that the Host Star has been found on the Star Catalogue
            transitHostStarMatch = true;
                                
            //REMINDER:
            //starCatalogue(starIterator, 0) = id(starIterator);
            //starCatalogue(starIterator, 1) = ra(starIterator);
            //starCatalogue(starIterator, 2) = dec(starIterator);
            //starCatalogue(starIterator, 3) = magn(starIterator);       
            
            
            //Set the ID position of the host source on the Star Catalogue on the DataSet
            hostStarID = starCatalogue(starIterator, 0);
            p_DataSet->datasetSethostStarID(hostStarID);

            //Set the magnitude of the host source on the DataSet
            hostStarMagnitude = starCatalogue(starIterator, 3);
            p_DataSet->datasetSethostStarMagnitude(hostStarMagnitude);
        

        } 
                  
        //Go for the next source in the Star catalogue
        starIterator++;
    }       
 
    //If there has NOT been any source in the Star Catalogue that matches with the Right Ascension and Declination of the transit planet 
    if (transitHostStarMatch == false) 
    {

        std::cerr << "\nError (ParamsTransit::ParamsTransitsearchHostStarOnCatalogue): Please, check that the Right Ascension and Declination of the transit planet corresponds with a source in the Star Catalogue with better accuracy than 0.00001 degrees." << std::endl;
        exit(1);
    }
           
}

//==============================================================================
