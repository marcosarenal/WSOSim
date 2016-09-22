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
        void ParamsTransitcalculation(DataSet &m_DataSet);
        void ParamsTransitsetting(DataSet &m_DataSet);
        void ParamsTransitcheckInputParameters(DataSet &m_DataSet);

private:
    
        DataSet*        p_DataSet;              //Pointer to the DataSet to retrieve parameters from it.

        Array<float, 1>  inTransitArray;        //Blitz array to be set into the DataSet.        
        
        bool   PerformExoTransit;               //Perform exoplanetary transit simulation for one source in the field (0=no/1=yes).
        double HostStarTransitRA;               //Right Ascension of the transit host star [deg] (It must match with a source in the star catalogue).
        double HostStarTransitDec;              //Declination of the transit host star [deg].
        double HostStarRadius;                  //Radius of the transit host star [Solar radius].
        double ExoplanetRadius;                 //Radius of the exoplanet [Solar radius].
        double ExoplanetOrbitalPeriod;          //Orbital period of the exoplanet [days].
        double PlanetaryOrbitSemiaxis;          //Semiaxis of the orbit of the exoplanet [AU].
        double PlanetaryOrbitInclination;       //Inclination of the orbit of the exoplanet as seen from Earth [deg].


};



#endif /* PARAMSTRANSIT_H */

