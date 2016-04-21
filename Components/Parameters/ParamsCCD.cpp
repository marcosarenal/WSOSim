///////////////////////////////////////////////////////////
//  ParamsCCD.cpp
//  Implementation of the Class ParamsCCD
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



#include "ParamsCCD.h"


//==============================================================================
/**
 * Constructor method
 */
ParamsCCD::ParamsCCD(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsCCD::~ParamsCCD(){}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================

/**
 * This function calls all the functions required to perform the CCD Preprocessing
 * calculation, including checking the validity of the input parameters, calculating the CCD
 * position and setting all the calculated parameters that are required in the
 * processing steps into the DataSet.
 * @param m_DataSet
 */
void ParamsCCD::paramsCCDcalculation(DataSet &m_DataSet)
{
    //Pointing to the dataSet
    p_DataSet = &m_DataSet;
    
    // Get the CCD parameters from the input parameters file.
    orientationFocalPlane  = p_DataSet->datasetGetorientationFocalPlane();
    opticalAxisDecCenter = p_DataSet->datasetGetopticalAxisDecCenter();
    opticalAxisRACenter = p_DataSet->datasetGetopticalAxisRACenter();
    
    
    //Resizing CCD for adding margins
    ParamsCCD::paramsCCDresizeCCD(m_DataSet);
    
    //Check whether the input parameters are valid
    ParamsCCD::paramsCCDcheckCCDParams(m_DataSet);
    ParamsCCD::paramsCCDmappingCCDParams();
    ParamsCCD::paramsCCDsubmappingCCDParams();
    
    ParamsCCD::paramsCCDPixelOffsetToRad();
    
    //Compute the orientation of the CCD without jitter
    ParamsCCD::paramsCCDcomputeCCDOrientation(0, 0, 0, 0, 0);
    ParamsCCD::paramsCCDcomputeRadiusFOVandCenter();
    
    //Set the modified parameters in DataSet
    p_DataSet->datasetSetCCDParams(ccdSizeX, ccdSizeY, subFieldSizeX, subFieldSizeY, originOffsetXmm,
                                   originOffsetYmm, raCenterSubField, declCenterSubField, radiusFOVCCD,
                                   edgePixels, raOpticalAxis, decOpticalAxis, rotationAngleOA, xFOVSubField, yFOVSubField);
    
    //Set the empty pixel map and subpixel map.
    p_DataSet->datasetSetinitPixelMap(initPixelMap);
    p_DataSet->datasetSetinitSubPixelMap(initSubPixelMap);
        
    p_DataSet->datasetSetpixelMap(pixelMap);
    p_DataSet->datasetSetsubPixelMap(subPixelMap);
}
//==============================================================================




//==============================================================================
/**
 * This function provides the predefined offsets for CCD.
 */
void ParamsCCD::paramsCCDresizeCCD(DataSet &m_DataSet)
{
    
    // Get the CCD parameters from the input parameters file.
    ccdSizeX = p_DataSet->datasetGetccdSizeX();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels();
    
    if (subPixelsPerPixel <= 0)
    {
            cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Number of SubPixels per Pixel must be > 0." << endl;
            exit(1);
    }
    halfPSFSize = (psfSubPixels/ subPixelsPerPixel + 1)/2;
    
    psfSize = 2 * halfPSFSize;
    
    ccdSizeX = ccdSizeX + psfSize;
    ccdSizeY = ccdSizeY + psfSize;
    
    // Get the Subfield size from the input parameters file.
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    
    subFieldSizeX = subFieldSizeX + psfSize;
    subFieldSizeY = subFieldSizeY + psfSize;

    
    
    LogManager::log<< "    Important: To provide proper treatment at the edges of the CCD and sub-image, the sizes of both are increased during the computations. ";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
    LogManager::log<< "    At the end of the computations the size of the sub-field is decreased again.";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
    LogManager::log<< "    At each edge, the CCD size is increased by " << halfPSFSize << " (in total " << psfSize
                   << " for each dimension)";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
    LogManager::log<< "    New CCD size (X, Y) = (" << ccdSizeX << ", " << ccdSizeY << ")";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
//	LogManager::log<< "New sub-field zero point (X, Y) = (" << subFieldZeroX << ", " << subFieldZeroY << ")";
//	GlobalVariables::logManager.LogManagerAppendLogAndShow();
    LogManager::log<< "    New sub-field size (X, Y) = (" << subFieldSizeX << ", " << subFieldSizeY << ")";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

    // The pointing direction of the CCD is not defined through its center of field but through the offset
    // of its origin (readout corner). The change in size of the CCD (above) to treat the edges correctly
    // must also result in a shift of the origin offset.
    // Here is computed the offset of the origin in order to get the CCD of normal size at the user defined position.
    // So, it is converted the halfPSFSize to mm; Then, we are in the focal plane, therefore everything is linear.
    
    // Get the Subfield size from the input parameters file.
    originOffsetXmm = p_DataSet->datasetGetOriginOffsetXmm();
    originOffsetYmm = p_DataSet->datasetGetOriginOffsetYmm();
    pixelSize = p_DataSet->datasetGetpixelSize();
    
    originOffsetXmm = originOffsetXmm - halfPSFSize * pixelSize / 1000.;
    originOffsetYmm = originOffsetYmm - halfPSFSize * pixelSize / 1000.;

    LogManager::log << "    New origin offset of enlarged CCD field (X, Y) = (" << originOffsetXmm << " mm, "
                    << originOffsetYmm << " mm)";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
 	
}
//==============================================================================
/**
 * This function checks that every CCD parameter is valid for performing the simulation.
 */
void ParamsCCD::paramsCCDcheckCCDParams(DataSet &m_DataSet)
{
    
	if (ccdSizeX <= 0 || ccdSizeY <= 0)
        
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Pixel dimensions of CCD must be > 0." << endl;
		exit(1);
	}
    
	if (subFieldSizeX <= 0 || subFieldSizeY <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Pixel dimensions of sub-field must be > 0." << endl;
		exit(1);
	}
    
    subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY = p_DataSet->datasetGetsubFieldZeroY();
	if (subFieldZeroX >= ccdSizeX || subFieldZeroY >= ccdSizeY || subFieldZeroX < 0 || subFieldZeroX < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Sub-field zeropoint must be inside CCD." << endl;
		exit(1);
	}
    
	if (subFieldZeroX + subFieldSizeX > ccdSizeX || subFieldZeroY + subFieldSizeY > ccdSizeY)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Sub-field pixel dimensions exceed size of CCD." << endl;
		exit(1);
	}
        
    
    areaTelescope = p_DataSet->datasetGetareaTelescope();
	if (areaTelescope <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Telescope Area must be > 0." << endl;
		exit(1);
	}
    
    transEff = p_DataSet->datasetGettransEff();
	if (transEff <= 0 || transEff > 1)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Transmission Efficiency must be 0 < t <= 1." << endl;
		exit(1);
	}
    
    quantEff = p_DataSet->datasetGetquantEff();
	if (quantEff <= 0 || quantEff > 1)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Quantum Efficiency must be 0 < qe <= 1." << endl;
		exit(1);
	}
    
    gain = p_DataSet->datasetGetgain();
	if (gain <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Gain must be > 0." << endl;
		exit(1);
	}
    
    pixelScale = p_DataSet->datasetGetpixelScale();
	if (pixelScale <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Pixel Scale must be > 0." << endl;
		exit(1);
	}
    
	if (pixelSize <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Pixel Size must be > 0." << endl;
		exit(1);
	}
    
	if (subPixelsPerPixel <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Number of SubPixels per Pixel must be > 0." << endl;
		exit(1);
	}
    
    fullWellSat = p_DataSet->datasetGetfullWellSat();
	if (fullWellSat <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Full-Well Saturation must be > 0." << endl;
		exit(1);
	}
    
    digitalSat = p_DataSet->datasetGetdigitalSat();
	if (digitalSat <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Digital Saturation must be > 0." << endl;
		exit(1);
	}
    
    fluxm0 = p_DataSet->datasetGetfluxm0();
	if (fluxm0 <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): m=0-Flux must be > 0." << endl;
		exit(1);
	}
    
	if (opticalAxisRACenter < 0 || opticalAxisRACenter > 360.)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Right Ascension of Optical Axis must be between 0 and 360." << endl;
		exit(1);
	}
    
	if (opticalAxisDecCenter < -90 || opticalAxisDecCenter > 90.)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Declination of Optical Axis must be between -90 and +90." << endl;
		exit(1);
	}
    
    readOutNoise = p_DataSet->datasetGetreadOutNoise();
	if (readOutNoise < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Readout Noise must be >= 0." << endl;
		exit(1);
	}
    
    readOutTime = p_DataSet->datasetGetreadOutTime();
	if (readOutTime < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Read-out Time must be >= 0." << endl;
		exit(1);
	}
    
    electronicOffset = p_DataSet->datasetGetelectronicOffset();
	if (electronicOffset < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Electronic Offset must be >= 0." << endl;
		exit(1);
	}
    
    numPrescanRows = p_DataSet->datasetGetnumPrescanRows();
	if (numPrescanRows < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Number of Prescan Rows must be >= 0." << endl;
		exit(1);
	}
    
    
    numSmearingOverscanRows = p_DataSet->datasetGetnumSmearingOverscanRows();
	if (numSmearingOverscanRows < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Number of Overscan Rows must be >= 0." << endl;
		exit(1);
	}
    
    ccdOrientation = p_DataSet->datasetGetccdOrientation();
    if (ccdOrientation < -360 || ccdOrientation > 360. )
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): CCD Orientation must be between -360 and 360." << endl;
		exit(1);
	}
    
    exposureTime = p_DataSet->datasetGetexposureTime();
	if (exposureTime <= 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Exposure Time must be > 0." << endl;
		exit(1);
	}
    
    flatfieldPixelNoise = p_DataSet->datasetGetflatfieldPixelNoise();
	if (flatfieldPixelNoise < 0 || flatfieldPixelNoise > 1)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Pixel-to-Pixel Noise must be 0 <= p <= 1." << endl;
		exit(1);
	}
    
    flatfieldWhiteNoise = p_DataSet->datasetGetflatfieldWhiteNoise();
	if (flatfieldWhiteNoise < 0 || flatfieldWhiteNoise > 1)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Sub-Pixel White Noise must be 0 <= p <= 1." << endl;
		exit(1);
	}
    
    
    meanCTE = p_DataSet->datasetGetmeanCTE();
	if (meanCTE < 0 || meanCTE > 1)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Mean CTE must be 0 <= CTE <= 1." << endl;
		exit(1);
	}
    
    cosmicHitRate = p_DataSet->datasetGetcosmicHitRate();
	if (cosmicHitRate < 0)
	{
		cerr << "\nError (ParamsCCD::paramsCCDcheckCCDParams): Cosmic Hit Rate must be >= 0." << endl;
		exit(1);
	}
    
    
}
//==============================================================================




//==============================================================================
/**
 *This function creates a initial pixel map with the subfield size.
 */
void ParamsCCD::paramsCCDmappingCCDParams()
{
    
    initPixelMap.resize(subFieldSizeX, subFieldSizeY);
    initPixelMap = 0.;	
    
    pixelMap.resize(subFieldSizeX, subFieldSizeY);
    pixelMap = 0.;
    
}
//==============================================================================




//==============================================================================
/**
 *The following function creates and fills with zeros a 2D matrix representing the lattice
 *of the CCD in subpixels i.e. each real pixel is subdivided in a sub-lattice of imaginary pixels.
 */
void ParamsCCD::paramsCCDsubmappingCCDParams()
{
    initSubPixelMap.resize(subPixelsPerPixel * subFieldSizeX, subPixelsPerPixel * subFieldSizeY);
    initSubPixelMap = 0.;	

    subPixelMap.resize(subPixelsPerPixel * subFieldSizeX, subPixelsPerPixel * subFieldSizeY);
    subPixelMap = 0.;
}
//==============================================================================




//==============================================================================
/**
 *This function converts the offset of the CCD from mm to radians with respect to the optical axis
 * this conversion does NOT (and SHOULD NOT) take into account the distortion of the projection of the sky to the focal plane since the
 * offset is measured only in the focal plane and we will compute the displacement in units of the focal plane
 * 	pixelScale ... arcsec/pixel
 * 	pixelSize ... size of pixel in micron (10^6 m)
 * 	originOffsetXmm ...  offset of the readout corner of the CCD in mm from the focal axis
 * 	ccdOriginOffsetX ... offset in radians
 * 	A=pixelScale
 * 	M=pixelSize
 *              A arcsec = 1 pixel
 *              1 pixel = M micron
 *              => A arcsec = M micron = 10^-6 m = 10^-3 mm
 *              1mm = 1000 micron = 1000 * A / M arcsec = 1000A/M * Pi/(3600*180) rad
 *              => N mm = N * 1000A/M *Pi/(3600/180) = N*3.6*A*Constants::Deg2Rad/M
 *
 */
void ParamsCCD::paramsCCDPixelOffsetToRad()
{
    ccdOriginOffsetX = originOffsetXmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);
    ccdOriginOffsetY = originOffsetYmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);
    
}
//==============================================================================




//==============================================================================
/**
 * General assumption: The telescope is pointed towards objects that have a infinite distance from the CCD. Therefore, any rotation of the satellite
 * has the same effect on a projection plane (CCD) independent of its distance to the jitter axis. Only the angular distance of a star from the
 * jitter axis is important.
 * This algorithm computes the position of a star on the CCD plane by taking into account the orientation of the optical axis, the jitter axis,
 * the jitter displacement (Euler) angles, and the position of the star on the sky. First, the jitter axes are defined such that z is aligned
 * with the axis of lowest inertia of the satellite and y is aligned such that the pitch results in N-S movement (y axis in equatorial plane).
 * Orientation of focal plane is defined by the user and oriented such that positive y-coordinates are towards N. The x-axis is in the equatorial plane.
 * The focal plane is rotated by jitter angles. Then, the star coordinates are transformed to carthesian coordinates and and projected gnomonic on the focal plane.
 * Finally, the position of the projected point relative to the CCD coordinate system is determined (through coordinate transformation).
 *      theta=0 and fi=0 is x,y,z=0,0,1
 * The rotation of the satellite (the jitter-axis-system) is defined as follows: A rotation about the yaw axis rotates also the other two axes. Next, the rotation
 * around the already rotated pitch-axis is carried out. This rotates the roll-axis. Finally, the rotation around the twice rotates roll-axis is performed.
 * Therefore, during each rotation of the CCD, the jitter-axes have to be rotated with it.
 *
 * IMPORTANT MODIFICATION on 2010 November 3:
 * The projection type has been changed from a normal projection to a GNOMONIC projection, since this type represents the actual projection
 * of an optical system much better. Two function in MathTools have been added that perform the gnomonic projection and its inverse.
 * The CCD is located somewhere (definable by the user) in the focal plane. The reference for the projection is the optical axis.
 *
 * determine orthogonal vectors defining the orientation of the Jitter axis (y=pitch toward N) in plane with given normal vector
 * roll around z
 * pitch around y
 * yaw around x
 *
 * @raJ @decJ  R.A. and Declination of jitter axis in degrees
 * @raC @decC  R.A. and Declination of CCD normal (=optical axis) in degrees
 * @yaw @pitch @roll  Jitter angles in arc sec
 */
void ParamsCCD::paramsCCDcomputeCCDOrientation(double raJitterAxis, double declJitterAxis, double yaw, double pitch, double roll)
{
    
	// Optical axis  orientation
	double Xcx, Xcy, Xcz; //vector in plane normal to ra,dec and in equatorial plane
	double Ycx, Ycy, Ycz; //vector in plane normal to y and z
	double Zcx, Zcy, Zcz; //normal vector of focal plane
    
    
	double fij = Constants::DEG2RAD * raJitterAxis;//jitter axis
	double thetaj = Constants::Pid2 - Constants::DEG2RAD * declJitterAxis; //convert to polar distance
	double fic = Constants::DEG2RAD * opticalAxisRACenter;//ccd orientation (normal)
	double thetac = Constants::Pid2 - Constants::DEG2RAD * opticalAxisDecCenter; //convert to polar distance
    
	//Orientation of jitter axes
	double Xjx, Xjy, Xjz; //vector in plane normal to y and z
	double Yjx, Yjy, Yjz; //vector in plane normal to ra,dec and in equatorial plane
	double Zjx, Zjy, Zjz; //normal vector of plane
    
	MathTools::polar2Karth(1., thetaj, fij, Zjx, Zjy, Zjz);
	MathTools::polar2Karth(1., Constants::Pid2, fij + Constants::Pid2, Yjx, Yjy, Yjz);
	//determine x axis
	MathTools::crossProduct(Yjx, Yjy, Yjz, Zjx, Zjy, Zjz, Xjx, Xjy, Xjz);

    
	//Orientation of optical (focal plane) axes
	MathTools::polar2Karth(1., thetac, fic, Zcx, Zcy, Zcz); //optical axis oriented towards the given coordinates
	MathTools::polar2Karth(1., Constants::Pid2, fic + Constants::Pid2, Ycx, Ycy, Ycz); //y-axis of focal plane. this can be oriented (rotated) as defined by the user
    
	//determine x axis = normal to y and z
	MathTools::crossProduct(Ycx, Ycy, Ycz, Zcx, Zcy, Zcz, Xcx, Xcy, Xcz);
    
	//Rotate field to match the user input for the orientation of the Focal Plane
	//Positive rotation is from N to W (same for jitter)
	//above, the focal plane is defined such that the y-axis points towards N
	MathTools::rotation(Xcx, Xcy, Xcz, Zcx, Zcy, Zcz, Constants::DEG2RAD * orientationFocalPlane);
	MathTools::rotation(Ycx, Ycy, Ycz, Zcx, Zcy, Zcz, Constants::DEG2RAD * orientationFocalPlane);
    

	//rotation angles are in arcsec. convert them to radians
	//yaw and pitch are inverted because the axes have also been inverted on the CCD with respect to the satellite axes
	yaw *= -Constants::DEG2RAD / 3600.;
	pitch *= -Constants::DEG2RAD / 3600.;
	roll *= Constants::DEG2RAD / 3600.;
    
	//now rotate the focal plane around the jitter axes
	//rotation is counterclockwise (positive pitch results generally in a negative delta(declination) on the CCD)
	//rotated focal plane coordinates in the order yaw, pitch and roll
	MathTools::rotation(Xcx, Xcy, Xcz, Xjx, Xjy, Xjz, yaw);
	MathTools::rotation(Ycx, Ycy, Ycz, Xjx, Xjy, Xjz, yaw);
	MathTools::rotation(Zcx, Zcy, Zcz, Xjx, Xjy, Xjz, yaw);
	//rotate also the other jitter axes
	MathTools::rotation(Zjx, Zjy, Zjz, Xjx, Xjy, Xjz, yaw);
	MathTools::rotation(Yjx, Yjy, Yjz, Xjx, Xjy, Xjz, yaw);
    
	//next rotation of focal plane around pitch axis
	MathTools::rotation(Xcx, Xcy, Xcz, Yjx, Yjy, Yjz, pitch);
	MathTools::rotation(Ycx, Ycy, Ycz, Yjx, Yjy, Yjz, pitch);
	MathTools::rotation(Zcx, Zcy, Zcz, Yjx, Yjy, Yjz, pitch);
	//rotate the jitter-roll axis
	MathTools::rotation(Zjx, Zjy, Zjz, Yjx, Yjy, Yjz, pitch);
    
	MathTools::rotation(Xcx, Xcy, Xcz, Zjx, Zjy, Zjz, roll);
	MathTools::rotation(Ycx, Ycy, Ycz, Zjx, Zjy, Zjz, roll);
	MathTools::rotation(Zcx, Zcy, Zcz, Zjx, Zjy, Zjz, roll);
    
	//Determine the orientation of the displaced focal plane
	double dummy;
	MathTools::karth2Polar(Zcx, Zcy, Zcz, dummy, decOpticalAxis, raOpticalAxis); //orientation of optical axis after displacement (jitter)
	//normal projection of the celestial pole (dec=90) into focal plane
	double px, py, pz; //projected pole in focal plane
	MathTools::normalProjectionPointOnPlane(0.f, 0.f, -1.f, Zcx, Zcy, Zcz, px, py, pz);
	//compute the angle between the displaced x-axis and the pole axis in the focal plane -> this should be the rotation of the field since the Gnomonic projection
	//has north always oriented at the top
	rotationAngleOA = MathTools::getAngle(px, py, pz, Xcx, Xcy, Xcz);
    
	//rotationAngleOA is between 0 and 180, so we must know the orientation of the two vectors
	//double nx, ny, nz;
	//MathTools::crossProduct(px, py, pz, Xcx, Xcy, Xcz, nx, ny, nz);
	//but I only require the x-component to determine if the orientation of the normal of p and X is in direction of Z or vice versa
	if ((py * Xcz - Xcy * pz) / Zcx < 0)
        {   
            rotationAngleOA = Constants::Pi2 - rotationAngleOA;
        }
    
}
//==============================================================================




//==============================================================================
/**
 * Computes the field of view in degrees from the dimension of the CCD and the pixel scale.
 */
void ParamsCCD::paramsCCDcomputeRadiusFOVandCenter()
{
	//each edge of the CCD has another angular dimension!
	//center of CCD: opticalAxisRACenter=>opticalAxisRACenter; opticalAxisDecCenter=>opticalAxisDecCenter
	ParamsCCD::paramsCCDgetTransformationCCDToSky(ccdSizeX / 2., ccdSizeY / 2., opticalAxisRACenter, opticalAxisDecCenter);
	//field of view of CCD
	double ra1, dec1, maxi, fov0;
	//int fieldExtension = 2 * Constants::DEG2RAD; //enlargement of the FOV to be sure that the complete CCD is covered with stars
    
	//Now determine the largest angular distance from the center of the CCD to each of its edges.
	//Each distance will be different due to the fact that the projection distorts the image
	ParamsCCD::paramsCCDgetTransformationCCDToSky(0, 0, ra1, dec1);
	maxi = MathTools::getAngularDistanceSphere(opticalAxisRACenter * Constants::DEG2RAD, opticalAxisDecCenter * Constants::DEG2RAD, ra1
                                               * Constants::DEG2RAD, dec1 * Constants::DEG2RAD);
	ParamsCCD::paramsCCDgetTransformationCCDToSky(ccdSizeX, 0, ra1, dec1);
	fov0 = MathTools::getAngularDistanceSphere(opticalAxisRACenter * Constants::DEG2RAD, opticalAxisDecCenter * Constants::DEG2RAD, ra1
                                               * Constants::DEG2RAD, dec1 * Constants::DEG2RAD);
	if (fov0 > maxi)
		maxi = fov0;
	ParamsCCD::paramsCCDgetTransformationCCDToSky(0, ccdSizeY, ra1, dec1);
	fov0 = MathTools::getAngularDistanceSphere(opticalAxisRACenter * Constants::DEG2RAD, opticalAxisDecCenter * Constants::DEG2RAD, ra1
                                               * Constants::DEG2RAD, dec1 * Constants::DEG2RAD);
	if (fov0 > maxi)
		maxi = fov0;
	ParamsCCD::paramsCCDgetTransformationCCDToSky(ccdSizeX, ccdSizeY, ra1, dec1);
	fov0 = MathTools::getAngularDistanceSphere(opticalAxisRACenter * Constants::DEG2RAD, opticalAxisDecCenter * Constants::DEG2RAD, ra1
                                               * Constants::DEG2RAD, dec1 * Constants::DEG2RAD);
	if (fov0 > maxi)
		maxi = fov0;
    
	radiusFOVCCD = Constants::RAD2DEG * maxi;
    
	//center of sub-field
	double x = (subFieldSizeX / 2. + subFieldZeroX);
	double y = (subFieldSizeY / 2. + subFieldZeroY);
    
	ParamsCCD::paramsCCDgetTransformationCCDToSky(x, y, raCenterSubField, declCenterSubField);
    
    //edgePixels Number of pixels by which the image size is increased on each side to ensure that the treatment of the edges is made correctly.
    //This should be half the size of the PSF.
    edgePixels = halfPSFSize;
    
	//fov of sub-field
	double ra2, dec2;
	ParamsCCD::paramsCCDgetTransformationCCDToSky(edgePixels + subFieldZeroX, edgePixels + subFieldZeroY + subFieldSizeY / 2., ra1, dec1);
	ParamsCCD::paramsCCDgetTransformationCCDToSky(edgePixels + subFieldZeroX + subFieldSizeX, edgePixels + subFieldZeroY + subFieldSizeY / 2., ra2,
                                                  dec2);
	xFOVSubField = Constants::RAD2DEG * MathTools::getAngularDistanceSphere(ra1 * Constants::DEG2RAD, dec1 * Constants::DEG2RAD, ra2
                                                                            * Constants::DEG2RAD, dec2 * Constants::DEG2RAD);
	ParamsCCD::paramsCCDgetTransformationCCDToSky(edgePixels + subFieldZeroX + subFieldSizeX / 2., edgePixels + subFieldZeroY, ra1, dec1);
	ParamsCCD::paramsCCDgetTransformationCCDToSky(edgePixels + subFieldZeroX + subFieldSizeX / 2., edgePixels + subFieldZeroY + subFieldSizeY, ra2,
                                                  dec2);
	yFOVSubField = Constants::RAD2DEG * MathTools::getAngularDistanceSphere(ra1 * Constants::DEG2RAD, dec1 * Constants::DEG2RAD, ra2
                                                                            * Constants::DEG2RAD, dec2 * Constants::DEG2RAD);
    
}
//==============================================================================




//==============================================================================
/**
 * Coordinates of star and gnomonic projection into the CCD focal plane.
 * @raS Star Right Ascension in degrees.
 * @decS Star Declination in degrees.
 */
void ParamsCCD::paramsCCDgetTransformationCCDToSky(double xStarCCD, double yStarCCD, double &raS, double &decS)
{
	double pixScalexRad = pixelScale * Constants::DEG2RAD / 3600.;
	xStarCCD *= pixScalexRad;
	yStarCCD *= pixScalexRad;
    
	double xRot, yRot;
	double shift = edgePixels * pixelScale * Constants::DEG2RAD / 3600.;
    
	MathTools::rotate(xStarCCD, yStarCCD, shift, shift, -ccdOrientation * Constants::DEG2RAD, xRot, yRot);
    
	MathTools::translate(xRot, yRot, -ccdOriginOffsetX, -ccdOriginOffsetY);
    
	MathTools::getInverseGnomonicProjection(xRot, yRot, raOpticalAxis, decOpticalAxis, -rotationAngleOA, raS, decS);
    
	raS *= Constants::RAD2DEG;
	decS *= Constants::RAD2DEG;
}
//==============================================================================
