///////////////////////////////////////////////////////////
//  StepTransit.h
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


#ifndef STEPTRANSIT_H
#define STEPTRANSIT_H

#include "DataSet.h"



/**
 * This class applies the modification in the flux of the source when the transit of the exoplanet takes place.
 * The reduction in flux is calculated as the fraction of the surface of the star ocultated by the surface of the planet.
 * The transit times are retrieved from the array inTransitArray with the beguining and end of each transit as in the example:
 * inTransitArray(initTransit1, endTransit1, initTransit2, endTransit2, ...) 
 */
class StepTransit
{

public:
	StepTransit();
	virtual ~StepTransit();
        void StepTransitapplication(DataSet &m_DataSet, double startTime);

        
private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
        double hostStarTransitRA;               //Right Ascension of the transit host star [deg] (It must match with a source in the star catalogue).
        double hostStarTransitDec;              //Declination of the transit host star [deg].
        double hostStarRadius;                  //Radius of the transit host star [Solar radius].
        double exoplanetRadius;                 //Radius of the exoplanet [Solar radius].

                
        Array<float, 1>  inTransitArray;        //Blitz array to be set into the DataSet.        
        Array<float, 2>  subPixelMap;              //Map array retrieved from DataSet.

};

#endif /* STEPTRANSIT_H */

