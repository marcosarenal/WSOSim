///////////////////////////////////////////////////////////
//  ParamsBackground.cpp
//  Implementation of the Class ParamsBackground
//  Created on:      05-Dic-2012 16:59:59 PM
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



#include "ParamsBackground.h"


//==============================================================================
/**
 * Constructor method
 */
ParamsBackground::ParamsBackground(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsBackground::~ParamsBackground(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * Sets the background flux at the center of the subfield of the CCD. The flux 
 * consists of zodiacal flux and diffuse stellar background (see class Sky).
 */
void ParamsBackground::ParamsBackgroundsetting(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    //Getting background from DataSet
    background = p_DataSet->datasetGetbackground();

    //if user set Background <0, compute and set the Background otherwise take the user value.
    if (background < 0) 
    {
        ParamsBackground::ParamsBackgroundcalculation(m_DataSet);
    }

    LogManager::log<< "    Sky background : "<< p_DataSet->datasetGetbackground()<< " e-/pixel/s ";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();	
}
//==============================================================================

//==============================================================================
/**
 * Compute the background flux at the center of the subfield of the CCD. The flux consists of zodiacal flux and diffuse stellar background (see class Sky).
 * The flux is in units of e- per normal pixel per second.
 * No image operation is made. To apply the flux to the image use the method apply(CCD ccd).
 * Sets the flux to -1 if the grid does not cover the position of the CCD.
 */
void ParamsBackground::ParamsBackgroundcalculation(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    Sky sky; // Sky class

    //TODO: Modify this wavelength dependence (pma)
    
    // Lower limit of wavelength interval     [m]
    lower_lambda = 400.0e-9;
    // Upper limit of wavelength interval     [m]
    upper_lambda = 800.0e-9; 

    ra = p_DataSet->datasetGetraCenterSubField() * Constants::DEG2RAD; 
    decl = p_DataSet->datasetGetdeclCenterSubField() * Constants::DEG2RAD; 

    flux0 = sky.ZodiacalFlux(lower_lambda, upper_lambda, ra, decl) + sky.StellarBgFlux(lower_lambda, upper_lambda, ra, decl);

    double lambda_c(630.0e-9);
    energy_photon = Constants::HPLANCK * Constants::CLIGHT / lambda_c;
    normalpixelarea = pow(Constants::DEG2RAD * p_DataSet->datasetGetpixelScale() / 3600., 2.);
    newBackground = (flux0 / energy_photon) * (p_DataSet->datasetGetareaTelescope() / 10000.) 
                * normalpixelarea * p_DataSet->datasetGettransEff() * p_DataSet->datasetGetquantEff();

    //Setting the background flux (in e-/pixel/s) calculated in the DataSet 
    p_DataSet->datasetSetbackground(newBackground); 

}
//==============================================================================
