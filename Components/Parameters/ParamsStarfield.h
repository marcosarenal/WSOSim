///////////////////////////////////////////////
//  ParamsStarfield.h
//  Implementation of the Class ParamsStarfield
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



#ifndef PARAMSSTARFIELD_H_
#define PARAMSSTARFIELD_H_


#include "DataSet.h"
#include "blitz/array.h"
#include "Constants.h"
#include "MathTools.h"
#include "FileUtilities.h"

#include <string>
#include <fstream>


using namespace blitz;

/**
 * Class which handles the star catalog that is used to create CCD/CMOS images. A
 * star catalog can be read in and used in a CCD object. A catalog must be an
 * ASCII file consisting of 3 columns separated by spaces: R.A, Dec and magnitude.
 */
class ParamsStarfield
{
    
public:
	ParamsStarfield();
	virtual ~ParamsStarfield();
    
    void ParamsStarfieldSetting(DataSet &m_DataSet);
    void ParamsStarfieldCoords(DataSet &m_DataSet);
    void ParamsStarfieldSetCCDParams(double a, double c, double d);
    void ParamsStarfieldReadCatalogue(string catalogueFileName);
    void ParamsStarfieldTransformationSkyToCCD(double raS, double decS, double &xStarCCD, double &yStarCCD);
    void ParamsStarfieldfillsubPixelMap(DataSet &m_DataSet);
    
    
    
private:
    
    
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.


    Array<float, 1> magn, ra, dec;
    Array<int,   1> id;


    Array<float, 2> subPixelStarListOnCCD, starListOnCCD;                   //Blitz map set into the DataSet.
    Array<float, 2> initSubPixelMap;                                        //Blitz map set into the DataSet.
    Array<float, 2> starCatalogue;                                          //Blitz map set into the DataSet.

    double radius_fov, ra_center, decl_center;
    double radiusFOVCCD, raCenterSubField, declCenterSubField;
    string catalogueFileName;                                               //Parameter retrieved from DataSet.
    int    NumStars, numStarsCCD, numStarsSubField;

    double pixelScale;                                                      //Parameter retrieved from DataSet.
    double originOffsetXmm, originOffsetYmm;                                //Parameter retrieved from DataSet.
    double originOffsetX, originOffsetY;                                    //Parameter retrieved from DataSet.
    int    edgePixels, pixelSize, subPixelsPerPixel;                        //Parameter retrieved from DataSet.
    double raOpticalAxis, decOpticalAxis, rotationAngleOA;                  //Parameter retrieved from DataSet.
    int    subFieldSizeX, subFieldSizeY;                                    //Parameter retrieved from DataSet.
    int    ccdSizeX, ccdSizeY;                                              //Parameter retrieved from DataSet.
    int    subFieldZeroX, subFieldZeroY;                                    //Parameter retrieved from DataSet.
    double ccdOrientation;                                                  //Parameter retrieved from DataSet.

    double exposureTime;                                                    //Parameter retrieved from DataSet.
    double fluxm0, areaTelescope, transEff, quantEff;                       //Parameter retrieved from DataSet.

    double xOpticalAxis, yOpticalAxis;                                      //Parameter retrieved from DataSet.
    
    
    
    
    
};
#endif /* PARAMSSTARFIELD_H_ */
