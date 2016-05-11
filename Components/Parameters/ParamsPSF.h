///////////////////////////////////////////////////////////
//  ParamsPSF.h
//  Implementation of the Class ParamsPSF
//  Created on:      23-Oct-2012 1:59:58 PM
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



#ifndef PARAMSPSF_H_
#define PARAMSPSF_H_

#include <string>
#include <fstream>

#include "DataSet.h"
#include "blitz/array.h"
#include "MathTools.h"
#include "Constants.h"
#include "FileUtilities.h"

//#include <opencv/cv.h>
//#include <opencv/cvwimage.h>

using namespace std;
using namespace blitz;


/**
 * This class performs the required computations for setting the PSF parameters in the DataSet.
 */
class ParamsPSF
{

public:

	ParamsPSF();
	virtual ~ParamsPSF();
    
    void paramsPSFCalculation(DataSet &m_DataSet);
    void paramsPSFcheckPSFParams();
    void paramsPSFgetPSFFileFromGrid(); 
    void paramsPSFComputeMask(DataSet &m_DataSet);
    void paramsPSFsetGaussian();
    void paramsPSFreadFromFile();
    void paramsPSFwriteSubPixelsToASCIIFile(string);
    void paramsPSFrotatePSF(Array<float, 2> &, double);

    
    static const int NO_PSF_ROTATION = 0;
	static const int PSF_ROTATION_OPTICAL_AXIS = 1;
	static const int PSF_ROTATION_ARBITRARY = 2;


private:

    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
    double angDist;                                         //Parameter retrieved from DataSet.
    double rotationAngle;                                   //Parameter retrieved from DataSet.
    double psfX;                                            //PSF orientation vector in x direction of pre-computed PSF.
    double psfY;                                            //PSF orientation vector in y direction of pre-computed PSF.
    int    psfRotation;                                     //Parameter retrieved from DataSet.
    double psfRotationAngle, psfOrientation;                //Parameter retrieved from DataSet.
    double raCenterSubField, declCenterSubField;            //Parameter retrieved from DataSet.
    double raOpticalAxis, decOpticalAxis;                   //Parameter retrieved from DataSet.
    double xOpticalAxis, yOpticalAxis;                      //Parameter retrieved from DataSet.
    
    string psfFileName;                                         //Parameter retrieved from DataSet.
    double psfGaussFWHM;                                    //Parameter retrieved from DataSet.
    int    subPixelsPerPixel;                               //Parameter retrieved from DataSet.
    bool   psfLocationDependent;                            //Parameter retrieved from DataSet.
    int    psfSubPixels, psfNumPixels;                      //Parameter retrieved from DataSet.
    double psfCenterX, psfCenterY, orientationAngle;        //Parameter retrieved from DataSet.
    int    psfSize;                                         //Parameter retrieved from DataSet.
    int    subFieldSizeX, subFieldSizeY;                    //Parameter retrieved from DataSet.
    int    subFieldZeroX, subFieldZeroY;                    //Parameter retrieved from DataSet.
    bool   useGauss;                                        //Parameter retrieved from DataSet.
    string outputPath;                                      //Parameter retrieved from DataSet.
    string prefix;                                          //Parameter retrieved from DataSet.
    string psfLocationFile;                                 //Parameter retrieved from DataSet.
    
    
    Array<float, 2> mask;                                   //Calculated PSF map to be set in the DataSet.

        
        
};
#endif /* PARAMSPSF_H_ */
