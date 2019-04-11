///////////////////////////////////////////////////////////
//  StepJitter.h
//  Implementation of the Class StepJitter
//  Created on:      23-Oct-2012 2:00:01 PM
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



#ifndef STEPJITTER_H_
#define STEPJITTER_H_


#include "LogManager.h"
#include "DataSet.h"
#include "FFT.h"
#include "Constants.h"
#include "MathTools.h"
#include "blitz/array.h"
#include "random/discrete-uniform.h"
#include "Statistics.h"

/**
 * This class applies given jitter Euler-angle displacements to a CCD object. A
 * jitter time series can be read from a file. A jitter time-series has to be
 * provided as an ASCII file which lists the pointing displacement in terms of the
 * Euler angles, yaw, pitch, and roll. To ensure a realistic modeling of jitter
 * variations, the time-step of the jitter time-series must be much smaller than
 * the integration time. The roll jitter axis of the space craft must not
 * necessarily be aligned with the optical axis of the telescope. The deviation
 * between these axes is fully taken into account for the calculation of the pixel
 * displacement of objects on the CCD. The definition of yaw, pitch, and roll is
 * such that if the jitter axis is aligned with the optical axis, the orientation
 * is 0 and the optical axis is pointed at Dec=0,  a positive yaw-angle results in
 * a displacement of the field toward East, a positive pitch results in a
 * displacement toward North, and a positive roll rotates the CCD clockwise (N-W-S-
 * E).
 */
class StepJitter
{
    
public:
	StepJitter();
	virtual ~StepJitter();
    void StepJitterapplication(DataSet &m_DataSet, double startTime);
    void StepJittercomputeFromFile(DataSet &m_DataSet, double startTime);
    void StepJittercomputeFromParameters(DataSet &m_DataSet, double startTime);
    void StepJittercomputeAxisOrientation();
    double StepJittercomputeCCDOrientation(double yaw, double pitch, double roll);
    void StepJitterCoords(double rotationAngleOA);
    void StepJitterTransformationSkyToCCD(double raS, double decS, double rotationAngleOA, double &xStarCCD, double &yStarCCD);
    void StepJitterfillsubPixelMap(double myExposure);
    
    
    
private:
    
    DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
    
    double raOpticalAxis;                                 //Parameter retrieved from DataSet. Orientation of optical axis after displacement (jitter) in radians
    double decOpticalAxis;                                //Parameter retrieved from DataSet. Orientation of optical axis after displacement (jitter) in radians
    double rotationAngleOA;                               //Parameter retrieved from DataSet. Rotation angle of optical axis after displacement (jitter) in radians
    double raJitterAxis, decJitterAxis;                   //Parameter retrieved from DataSet.
    bool   useJitterFromFile;                             //Parameter retrieved from DataSet.
    double jitterPosAngle;                                //Parameter retrieved from DataSet.
    double jitterMultFactor;                              //Parameter retrieved from DataSet.
    double jitterAngularDist;                             //Parameter retrieved from DataSet.
    double exposureTime, readOutTime;                     //Parameter retrieved from DataSet.
    double integrationTime;                               //integrationTime = exposureTime + readOutTime;
    
    int    NumStars, numStarsCCD, numStarsSubField;
    
    double pixelScale;                                    //Parameter retrieved from DataSet.
    double originOffsetXmm, originOffsetYmm;              //Parameter retrieved from DataSet.
    double originOffsetX, originOffsetY;                  //Parameter retrieved from DataSet.
    int    edgePixels, pixelSize, subPixelsPerPixel;      //Parameter retrieved from DataSet.
    int    ccdSizeX, ccdSizeY, subFieldZeroX;             //Parameter retrieved from DataSet.
    int    subFieldZeroY, subFieldSizeX, subFieldSizeY;   //Parameter retrieved from DataSet.
    double ccdOrientation;                                //Parameter retrieved from DataSet.
    double fluxm0, areaTelescope, transEff, quantEff;     //Parameter retrieved from DataSet.
	double opticalAxisRACenter, opticalAxisDecCenter;     //Parameter retrieved from DataSet.
	double orientationFocalPlane;                         //Parameter retrieved from DataSet.
     
    double jitterInterval;                                //Parameter retrieved from DataSet.
    double jitterRms;                                     //Parameter retrieved from DataSet.
    double jitterDrift;                                   //Parameter retrieved from DataSet.
    double jitterRepointing;                              //Parameter retrieved from DataSet.
    double yawJitter;                                     //Random generated parameter.
    double pitchJitter;                                   //Random generated parameter.
    double rollJitter;                                    //Random generated parameter.
        
    double yawPrevious;
    double pitchPrevious;
    double rollPrevious; 
            
    Array<double, 2> jitterInputParams;                   //Map array retrieved from DataSet.
    Array<float, 2>  subPixelMap;                         //Map array retrieved from DataSet.
    Array<float, 2>  starListOnCCD;                       //Map array retrieved from DataSet.
    Array<float, 2>  jitterMap;                           //Map array updated into DataSet.
    Array<float, 2>  subPixelStarListOnCCD;               //Map array retrieved from DataSet.
    Array<float, 2>  starCatalogue;                       //Array retrieved from DataSet.

};
#endif /* STEPJITTER_H_ */