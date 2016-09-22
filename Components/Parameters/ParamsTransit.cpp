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
    HostStarTransitRA = p_DataSet->datasetGethostStarTransitRA();
    HostStarTransitDec = p_DataSet->datasetGethostStarTransitDec();
    HostStarRadius = p_DataSet->datasetGethostStarRadius();
    ExoplanetRadius = p_DataSet->datasetGetexoplanetRadius();
    ExoplanetOrbitalPeriod = p_DataSet->datasetGetexoplanetOrbitalPeriod();
    PlanetaryOrbitSemiaxis = p_DataSet->datasetGetplanetaryOrbitSemiaxis();
    PlanetaryOrbitInclination = p_DataSet->datasetGetplanetaryOrbitInclination();
        
    //Check whether the input parameters are valid
    ParamsTransit::ParamsTransitcheckInputParameters(m_DataSet);

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
    PerformExoTransit = p_DataSet->datasetGetperformExoTransit();
	if (PerformExoTransit)
	{

            if (HostStarTransitRA < 0 || HostStarTransitRA > 360)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Right Ascension of the transit host star must be between 0 and 360." << endl;
		exit(1);
            }

            if (HostStarTransitDec < -90 || HostStarTransitDec > 90)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Declination of the transit host star must be between -90 and +90." << endl;
		exit(1);
            }
            
            if (HostStarRadius < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the host star must be > 0." << endl;
		exit(1);
            }
            
            if (ExoplanetRadius < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Radius of the transit planet  must be > 0." << endl;
		exit(1);
            }
                        
            if (HostStarRadius < ExoplanetRadius )
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Please, check that the radius of the transit planet is smaller than the radius of its host star." << endl;
		exit(1);
            }

            if (ExoplanetOrbitalPeriod < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Exoplanet Orbital Period must be > 0." << endl;
		exit(1);
            }
            
            if (PlanetaryOrbitSemiaxis < 0)
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Semiaxis must be > 0." << endl;
		exit(1);
            }

            if (PlanetaryOrbitInclination < 0 || PlanetaryOrbitInclination > 180 )
            {
		cerr << "\nError (ParamsTransit::ParamsTransitcheckInputParameters): Planetary Orbit Inclination must be between 0 and 180." << endl;
		exit(1);
            }


	}
    
    
}
//==============================================================================

//==============================================================================
/**
 * Compute the background flux at the center of the subfield of the CCD. The flux consists of zodiacal flux and diffuse stellar background (see class Sky).
 * The flux is in units of e- per normal pixel per second.
 * No image operation is made. To apply the flux to the image use the method apply(CCD ccd).
 * Sets the flux to -1 if the grid does not cover the position of the CCD.
 */
void ParamsTransit::ParamsTransitcalculation(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;



}
//==============================================================================
