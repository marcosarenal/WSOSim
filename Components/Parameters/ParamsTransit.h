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




#ifndef PARAMSTRANSIT_H
#define PARAMSTRANSIT_H

#include "GlobalVariables.h"
#include "DataSet.h"
#include "Constants.h"
#include "MathTools.h"
#include <math.h>

/**
 * This class computes and sets the time where the exoplanet transit takes place. These values
 * are set in the array inTransitArray as in the example:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...)
 * 
 */
class ParamsTransit
{

public:
	ParamsTransit();
	virtual ~ParamsTransit();
        void ParamsTransittimeCalculation(DataSet &m_DataSet);
        void ParamsTransitsetting(DataSet &m_DataSet);
        void ParamsTransitcheckInputParameters(DataSet &m_DataSet);
        void ParamsTransitsearchHostStarOnCatalogue(DataSet &m_DataSet);
        
private:
    
        DataSet*        p_DataSet;              //Pointer to the DataSet to retrieve parameters from it.

        Array<float, 1>  inTransitArray;        //Blitz array to be set into the DataSet.        
        Array<float, 2> starCatalogue;          //Blitz map set into the DataSet.

        
        bool   performExoTransit;               //Perform exoplanetary transit simulation for one source in the field (0=no/1=yes).
        double hostStarTransitRA;               //Right Ascension of the transit host star [deg] (It must match with a source in the star catalogue).
        double hostStarTransitDec;              //Declination of the transit host star [deg].
        double hostStarRadius;                  //Radius of the transit host star [Solar radius].
        double exoplanetRadius;                 //Radius of the exoplanet [Solar radius].
        double exoplanetOrbitalPeriod;          //Orbital period of the exoplanet [days].
        double planetaryOrbitSemiaxis;          //Semiaxis of the orbit of the exoplanet [AU].
        double planetaryOrbitInclination;       //Inclination of the orbit of the exoplanet as seen from Earth [deg].

        int    numExposures;                   //Parameter retrieved from DataSet.
        double exposureTime;                   //Parameter retrieved from DataSet.
        double readOutTime;                    //Parameter retrieved from DataSet.
        double integrationTime;                //integrationTime = exposureTime + readOutTime.
        double startTime;                      //Starting time for each exposure to be applied the Jitter
        double totalSimulationDuration;        //Duration of the whole simulation [in seconds]
        double orbitalPeriodInSeconds;         //Orbital period [in seconds]
        int    numberOfTransits;               //Number of transits in the whole simulation.
        double impactParameter;                //impact parameter is the sky-projected distance between the centre of the stellar disc and the centre of the planetary disc at conjunction
        double transitDuration;                //Transit duration [s]
        double planetaryOrbitSemiaxisInSolarRadius;                //planetary Orbit Semiaxis In Solar Radius units
        int    hostStarID;                     //ID position of the host source on the Star Catalogue
        double hostStarMagnitude;              //Magnitude of the host source on the Star Catalogue
        
};



#endif /* PARAMSTRANSIT_H */

