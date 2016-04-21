///////////////////////////////////////////////////////////
//  ParamsBackground.h
//  Implementation of the Class ParamsBackground
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



#ifndef PARAMSBACKGROUND_H_
#define PARAMSBACKGROUND_H_

#include "GlobalVariables.h"
#include "DataSet.h"
#include "Sky.h"


/**
 * This class computes and sets the sky background value on the detector for a certain 
 * location in the sky creating a 2D array in the pixel coordinates that will take into 
 * account the additional flux from the zodiacal light, the galactic and extra galactic 
 * light and the unresolved stars light.
 */
class ParamsBackground
{

public:
	ParamsBackground();
	virtual ~ParamsBackground();
        void ParamsBackgroundcalculation(DataSet &m_DataSet);
        void ParamsBackgroundsetting(DataSet &m_DataSet);

private:
    
        DataSet*        p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.

    	double lower_lambda, upper_lambda;   
	double ra, decl;
	double flux0;                            
        double energy_photon;
        double normalpixelarea;
        double background;                      //Parameter retrieved from DataSet.
        float  newBackground;                   //Parameter set into DataSet.



};
#endif /* PARAMSBACKGROUND_H_ */

