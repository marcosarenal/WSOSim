///////////////////////////////////////////////////////////
//  StepJitter.cpp
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



#include "StepJitter.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepJitter::StepJitter(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepJitter::~StepJitter(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function applies the Jitter to the imaging. 
 * The jitter to be applied might be calculated using an algorithm and input parameters or just 
 * read from an input file.
 * @param m_DataSet Reference to the object DataSet containing all the simulation parameters
 * @param startTime Start time of the exposure to be applied the jitter
 */
void StepJitter::StepJitterapplication(DataSet &m_DataSet, double startTime)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    raOpticalAxis = p_DataSet->datasetGetraOpticalAxis();       
    decOpticalAxis = p_DataSet->datasetGetdecOpticalAxis();     
    rotationAngleOA = p_DataSet->datasetGetrotationAngleOA();   
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    
    useJitterFromFile = p_DataSet->datasetGetuseJitterFromFile();
    jitterAngularDist = p_DataSet->datasetGetjitterAngularDist();
    jitterPosAngle = p_DataSet->datasetGetjitterPosAngle();
    jitterMultFactor = p_DataSet->datasetGetjitterMultFactor();
    
    //Getting exposureTime an readOutTime (in secs) from DataSet
    exposureTime = p_DataSet->datasetGetexposureTime();
    readOutTime = p_DataSet->datasetGetreadOutTime();
    integrationTime = exposureTime + readOutTime;

    //Retrieving parameters from DataSet
    fluxm0 = p_DataSet->datasetGetfluxm0();
    areaTelescope  = p_DataSet->datasetGetareaTelescope();
    transEff = p_DataSet->datasetGettransEff();
    quantEff = p_DataSet->datasetGetquantEff();
    ccdSizeX = p_DataSet->datasetGetccdSizeX();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();
    subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY = p_DataSet->datasetGetsubFieldZeroY();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    ccdOrientation = p_DataSet->datasetGetccdOrientation();
            
    //Getting detector parameters
    pixelScale = p_DataSet->datasetGetpixelScale();
    originOffsetXmm = p_DataSet->datasetGetOriginOffsetXmm();
    originOffsetYmm = p_DataSet->datasetGetOriginOffsetYmm();
    pixelSize = p_DataSet->datasetGetpixelSize();
       
    originOffsetX = originOffsetXmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);
    originOffsetY = originOffsetYmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);
        
    orientationFocalPlane  = p_DataSet->datasetGetorientationFocalPlane();
    opticalAxisDecCenter = p_DataSet->datasetGetopticalAxisDecCenter();
    opticalAxisRACenter = p_DataSet->datasetGetopticalAxisRACenter();
    

    //The subPixelMap is initialized for each exposure. 
    //Initializing the subpixel map from DataSet (it is null at this point)
    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
    subPixelMap = 0.0;
    
    // Compute the orientation (Ra, Dec) of an axis from the given CCD center of field
    // coordinates (RA, Dec) and the angular distance of the axis and its position angle.
    StepJitter::StepJittercomputeAxisOrientation();
                

    //If the jitter pointing positions are provided from an input file
    if (useJitterFromFile)
    {    
        
        //Compute JITTER using the jitter pointing input file
        StepJitter::StepJittercomputeFromFile(m_DataSet, startTime);        
    }
    
    //If the jitter pointing positions are not provided from an input file... compute from parameters
    else
    {
        //... this function calculate the jitter pointing positions using input parameters
        StepJitter::StepJittercomputeFromParameters(m_DataSet, startTime);
    }
    
                
    
    //Setting the new calculated subPixelMap into the DataSet
    m_DataSet.datasetSetsubPixelMap(subPixelMap);  

    
    LogManager::log <<"    Successfully added Jitter.";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================





//==============================================================================
/**
 * This function iterates in the different jitter position provided from a jitter file. 
 * Here is computed the jitter image by computing the star positions from the starlist
 * and shifting them. The setting in the subPixelMap is performed again but, in
 * this case, for each of the jittering positions in the jitter Input positions file.
 * The convolution is made afterwards in the StepConvolvePSFapplication module.
 * @param startTime Start time of the exposure to be applied the jitter

 */
void StepJitter::StepJittercomputeFromFile(DataSet &m_DataSet, double startTime)
{
            
    //Retrieving the jitter Input Parameters from DataSet
    //Initialize jitterInputParams array  
    jitterInputParams.resize(m_DataSet.datasetGetJitterInputParams().extent(0),m_DataSet.datasetGetJitterInputParams().extent(1));
    jitterInputParams = 0.0;
    jitterInputParams = m_DataSet.datasetGetJitterInputParams();
    
    //Return when the startTime reach the end of the time-serie in the jitter file
    if (startTime > max(jitterInputParams(Range::all(), 0)))
    {

		std::cerr << "\nError (StepJitter::StepJittercomputeFromFile()): Reached the end of the time-serie in the jitter file before the end of the simulated time-serie." << std::endl;
		std::cerr << "Try a larger jitter input file, a smaller number of simulations, or use the jitter parameters." << std::endl;
        exit(1);
    }
      

	//Initialize time and time interval between jitter positions
    double dt = 0.0;
    double  t = 0.0;
                
    int iterJitterPointing = 0;
    
    //Increase the counter in the jitterInputParams to match the exposure startTime 
    //with the time in the jitterInputParams
    while (jitterInputParams(iterJitterPointing, 0) < startTime)
	{	 
		iterJitterPointing++;
	}

	//When startTime overtakes the time in the jitterInputParams, back one unit  
	if (startTime < jitterInputParams(iterJitterPointing, 0))
    {
		iterJitterPointing--;
    }
	                
    //Loop during the exposure time adding jitter intervals time (dt = jitterInputParams(iterJitterPointing + 1, 0) - jitterInputParams(iterJitterPointing, 0)) 
	do
    {
        //If jitter interval is larger than exposure time
        if ( exposureTime < jitterInputParams(iterJitterPointing + 1, 0) - jitterInputParams(iterJitterPointing, 0))
        {
            dt = exposureTime;            

        }
        //To integrate in time from the beginning of the exposure to the first jitter position (if smaller than exposure time)
        else if (t == 0 && startTime > jitterInputParams(iterJitterPointing, 0))
        {
            dt = jitterInputParams(iterJitterPointing + 1, 0) - startTime;            

        }
        //To integrate in time from one jitter position to the next one (if smaller than exposure time)
        else if (t + exposureTime <= exposureTime)
        {
            dt = jitterInputParams(iterJitterPointing + 1, 0) - jitterInputParams(iterJitterPointing, 0);
        }
        //to integrate in time from the last jitter position to the end of this exposure 
        else
        {
            dt = exposureTime - t;
        }

        // yaw, pitch and roll are given as input to compute the CCD orientation for each jitter position
        rotationAngleOA = StepJitter::StepJittercomputeCCDOrientation(jitterInputParams(iterJitterPointing, 1), jitterInputParams(iterJitterPointing, 2), jitterInputParams(iterJitterPointing, 3));
		             

         //Each star position in the detector is recalculated according to the jitter position and overwritten in starListOnCCD
        StepJitter::StepJitterCoords(rotationAngleOA);
        		 

        //Each jitter is summed in the subPixelMap
        StepJitter::StepJitterfillsubPixelMap(dt);     
                
            
		iterJitterPointing++;
		t = t + dt;

	} while (t < exposureTime && iterJitterPointing < jitterInputParams.rows() - 1);
    
    		 
        
}
//==============================================================================
    




//==============================================================================
/**
 * This function iterates in the different jitter position generating the jitter 
 * positions apart from a generated distribution. 
 * @param startTime Start time of the exposure to be applied the jitter.
 */
void StepJitter::StepJittercomputeFromParameters(DataSet &m_DataSet, double startTime)
{
    //input parameters:
    jitterInterval = p_DataSet->datasetGetjitterInterval(); // in seconds
    jitterRms = p_DataSet->datasetGetjitterRms();      // in arcsecs
    jitterDrift = p_DataSet->datasetGetjitterDrift(); //  in arcsec/min
        
    //Every jitterRepointing time [hours], point back to the initial position (center of the FoV)
    jitterRepointing = p_DataSet->datasetGetjitterRepointing();  

    double repointingCadence = jitterRepointing * 3600 ; //from hours to seconds
    
    //Jitter pointing in yaw, pitch and roll described through normally distributed random variable
    // with zero mean and jitterRms standard deviation
    Normal<double, ranlib::MersenneTwister, ranlib::independentState> yawRand(0, jitterRms);
    Normal<double, ranlib::MersenneTwister, ranlib::independentState> pitchRand(0, jitterRms);
    Normal<double, ranlib::MersenneTwister, ranlib::independentState> rollRand(0, jitterRms);

    //Seed differently these jitter pointings.
    yawRand.seed(DataSet::seedRNG + startTime + 34);
    pitchRand.seed(DataSet::seedRNG+ startTime + 93);
    rollRand.seed(DataSet::seedRNG+ startTime + 73); 
        
       

	//Initialize time and time interval between jitter positions
    double t = 0.0;
             
   
    //Get in here to retrieve the previous pointing except for the first time (startTime = 0)
    //or after the re-pointing
    if (fmod(startTime + integrationTime,repointingCadence) > fmod(integrationTime,repointingCadence)   && (startTime !=0))
    {    
        //Loop in the exposure time to set (jitterInterval / exposureTime) times the flux in the subPixelMap
        while(t < exposureTime)
        {
        
            //Each tilt moves apart from the previous position + Jitter  + drift
            yawJitter   =  yawPrevious   + yawRand.random() + (jitterDrift/60) * t ;
            pitchJitter =  pitchPrevious + pitchRand.random() + (jitterDrift/60) * t;
            rollJitter  =  rollPrevious  + rollRand.random()  + (jitterDrift/60) * t;


            // yaw, pitch and roll are given as input to compute the CCD orientation for each jitter position
            rotationAngleOA = StepJitter::StepJittercomputeCCDOrientation(yawJitter, pitchJitter, rollJitter);


            //Each star position in the detector is recalculated according to the jitter position and overwritten in starListOnCCD
            StepJitter::StepJitterCoords(rotationAngleOA);


            //Each jitter is summed in the subPixelMap
            StepJitter::StepJitterfillsubPixelMap(jitterInterval);  
            t = t + jitterInterval;
        }
    }
    
    //First jitter (startTime = 0) or after re-pointing to the initial position.
    //Here is given the initial position and calculated the jitter
    else
    {
        
        while(t < exposureTime)
        {
            //Point to origin in the beginning of the observation or after repointing
            if (t == 0)
            {
                //yaw, pitch and roll initial positions
                yawPrevious =  0.0;
                pitchPrevious = 0.0;
                rollPrevious = 0.0 ;

                yawJitter   = 0.0;
                pitchJitter = 0.0;
                rollJitter  = 0.0;
                
            }

                       
            // yaw, pitch and roll are given as input to compute the CCD orientation for each jitter position
            rotationAngleOA = StepJitter::StepJittercomputeCCDOrientation(yawJitter, pitchJitter, rollJitter);

            //Each star position in the detector is recalculated according to the jitter position and overwritten in starListOnCCD
            StepJitter::StepJitterCoords(rotationAngleOA);

            //Each jitter is summed in the subPixelMap
            StepJitter::StepJitterfillsubPixelMap(jitterInterval);  
                   
            t = t + jitterInterval;
            
        }
        
        
    }
    
    //Set the current yaw, pitch and roll values for the next exposure
    yawPrevious   = yawPrevious  + (jitterDrift/60) * exposureTime;
    pitchPrevious = pitchPrevious  + (jitterDrift/60) * exposureTime;
    rollPrevious  = rollPrevious  + (jitterDrift/60) * exposureTime ;
     
 
//    LogManager::log <<"    Successfully added Jitter from parameters.";
//    GlobalVariables::logManager.LogManagerShowLog();
}
//==============================================================================





//==============================================================================
/**
 * Compute the orientation (Ra, Dec) of an axis from the given CCD center of field
 * coordinates (RA, Dec) and the angular distance of the axis and its position angle.
 * All parameters in degrees.
 * Problem still to solve: if one of the axes points towards the north pole, computations are not correct anymore
 */
void StepJitter::StepJittercomputeAxisOrientation()
{
    
	double fic =  raOpticalAxis;
	double thetac = Constants::Pid2 -  decOpticalAxis; //convert to polar distance
    
    
	//Orientation of optical axis
	double Xcx, Xcy, Xcz; //vector in plane normal to ra,dec and in equatorial plane
	double Zcx, Zcy, Zcz; //normal vector of CCD plane
	MathTools::polar2Karth(1.f, thetac, fic, Zcx, Zcy, Zcz);
	MathTools::polar2Karth(1.f, Constants::Pid2, fic + Constants::Pid2, Xcx, Xcy, Xcz);
    
	//rotate around the optical axis with the given jitter axis orientation angles to find the true jitter axis
	double Zjx = Zcx;
	double Zjy = Zcy;
	double Zjz = Zcz;
	MathTools::rotation(Zjx, Zjy, Zjz, Xcx, Xcy, Xcz, -jitterAngularDist * Constants::DEG2RAD); //rotation around the x-axis
	MathTools::rotation(Zjx, Zjy, Zjz, Zcx, Zcy, Zcz, -jitterPosAngle * Constants::DEG2RAD); //rotation around the z-axis
    
    
	//transform jitter axis coordinates to spherical system
	double radius, thetaj, fij;
	MathTools::karth2Polar(Zjx, Zjy, Zjz, radius, thetaj, fij);
    
	raJitterAxis = fij * Constants::RAD2DEG;
	decJitterAxis = thetaj * Constants::RAD2DEG;
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
 * @param raJ decJ  R.A. and Declination of jitter axis in degrees
 * @param raC decC  R.A. and Declination of CCD normal (=optical axis) in degrees
 * @param yaw pitch roll  Jitter angles in arc sec
 * @return rotationAngleOA Rotation angle of the optical axis for each jitter position (in radians).
 */
double StepJitter::StepJittercomputeCCDOrientation(double yaw, double pitch, double roll)
{
    
	// Optical axis  orientation
	double Xcx, Xcy, Xcz; //vector in plane normal to ra,dec and in equatorial plane
	double Ycx, Ycy, Ycz; //vector in plane normal to y and z
	double Zcx, Zcy, Zcz; //normal vector of focal plane
    
    
	double fij = Constants::DEG2RAD * raJitterAxis;//jitter axis
	double thetaj = Constants::Pid2 - Constants::DEG2RAD * decJitterAxis; //convert to polar distance
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
    
    if (rotationAngleOA > Constants::Pi2)
    {

        rotationAngleOA =  rotationAngleOA - Constants::Pi2;
    }
    
//    if (rotationAngleOA > Constants::Pi)
//    {
//
//        rotationAngleOA =  rotationAngleOA - Constants::Pi;
//    }
    
   
    return rotationAngleOA;
    
}
//==============================================================================




//==============================================================================
/**
 * Compute the x,y-coordinates on the detector of the stars in the starfield.
 * This calculation has already be made in the ParamsStarfield module but when the 
 * jittering is active the calculation has to be performed for each of the jitter positions.
 * 
 * @param rotationAngleOA Rotation angle of the optical axis for each jitter position (in radians).
 */
void StepJitter::StepJitterCoords(double rotationAngleOA)
{
    // Retrieving the starCatalogue: contains for each star, its X and Y position in CCD, magnitude, RA,
    // declination and identification number.
    starCatalogue.resize(p_DataSet->datasetGetStarCatalogue().extent(0),p_DataSet->datasetGetStarCatalogue().extent(1));
    starCatalogue = 0.0;
    starCatalogue = p_DataSet->datasetGetStarCatalogue();
 
    
    //Initialize for each jitter iteration
    NumStars = starCatalogue.rows();
    numStarsCCD = 0;
    numStarsSubField = 0;
     
    //starListonCCD and subPixelStarListOnCCD will be overwritten in the DataSet 
    starListOnCCD.resize(NumStars, 6);
    starListOnCCD = 0.0;
    
    subPixelStarListOnCCD.resize(NumStars, 4);
    subPixelStarListOnCCD = 0.0;
    
    double xTrue = 0, yTrue = 0;
    int xTrueInt, yTrueInt;
            
    //Iteration in each star in the starCatalogue to check whether it is on the CCD and calculate its position if it is.
    for (int i = 0; i < NumStars; i++)
    { 
        //This function takes right ascension and declination of each star in the starCatalogue and the rotation angle
        //of the optical axis to get its position in the detector.
        StepJitter::StepJitterTransformationSkyToCCD(starCatalogue(i, 1), starCatalogue(i, 2), rotationAngleOA, xTrue, yTrue);
                                    
        //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at
        //intrapixel level (because we round and int)
        xTrueInt = int(round(xTrue * subPixelsPerPixel));
        yTrueInt = int(round(yTrue * subPixelsPerPixel));

       
        //For each jitter position there is a new position of each star in the detector. Here will be overwritten
        if (xTrueInt >= 0 && yTrueInt >= 0 && xTrueInt < ccdSizeX * subPixelsPerPixel && yTrueInt < ccdSizeY * subPixelsPerPixel)
        {
            //star is in the FoV of the CCD
            starListOnCCD(numStarsCCD, 0) = xTrue * subPixelsPerPixel; //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at
            starListOnCCD(numStarsCCD, 1) = yTrue * subPixelsPerPixel; //intrapixel level (because we round and int)
            starListOnCCD(numStarsCCD, 2) = starCatalogue(i, 3);  //magnitude
            starListOnCCD(numStarsCCD, 3) = starCatalogue(i, 1);  //RA
            starListOnCCD(numStarsCCD, 4) = starCatalogue(i, 2);  //declination
            starListOnCCD(numStarsCCD, 5) = starCatalogue(i, 0);  //star ID

            numStarsCCD++;            
            
            if (xTrueInt >= subFieldZeroX * subPixelsPerPixel && yTrueInt >= subFieldZeroY * subPixelsPerPixel && xTrueInt
                < (subFieldZeroX + subFieldSizeX) * subPixelsPerPixel && yTrueInt < (subFieldZeroY + subFieldSizeY)
                * subPixelsPerPixel)
            {
                //star is in the FoV of the sub image
                subPixelStarListOnCCD(numStarsSubField, 0) = (xTrue - subFieldZeroX) * subPixelsPerPixel; //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at
                subPixelStarListOnCCD(numStarsSubField, 1) = (yTrue - subFieldZeroY) * subPixelsPerPixel; //intrapixel level (because we round and int)
                subPixelStarListOnCCD(numStarsSubField, 2) = starCatalogue(i, 3);//magnitude
                subPixelStarListOnCCD(numStarsSubField, 3) = starCatalogue(i, 0) ;//star ID
                    
                numStarsSubField++;
     
            }
        }    
    }
    
    
    // Defining the starListOnCCD blitz array
    Array<double, 2> subPixelStarListOnCCDResized;
    Array<double, 2> starListOnCCDResized;

    //Remove zeros from the starListOnCCD array
    starListOnCCDResized.resize(numStarsCCD, 6);
    starListOnCCDResized = 0.0;
    starListOnCCDResized(Range::all(),Range::all()) = starListOnCCD(Range(0,numStarsCCD), Range::all());

    starListOnCCD.resize(numStarsCCD, 6);
    starListOnCCD = 0.0;
    starListOnCCD = starListOnCCDResized;


    //Remove zeros from the subPixelStarListOnCCD array
    subPixelStarListOnCCDResized.resize(numStarsSubField, 4);
    subPixelStarListOnCCDResized = 0.0;
    subPixelStarListOnCCDResized(Range::all(), Range::all()) = subPixelStarListOnCCD(Range(0,numStarsSubField), Range::all());


    subPixelStarListOnCCD.resize(numStarsSubField, 4);
    subPixelStarListOnCCD = 0.0;
    subPixelStarListOnCCD = subPixelStarListOnCCDResized;
        
    
    //Setting the star list in the DataSet
    p_DataSet->datasetSetStarListOnCCD(starListOnCCD);
    
    //Setting the star list in the DataSet to be used in StepJitter
    p_DataSet->datasetSetsubPixelStarListOnCCD(subPixelStarListOnCCD);
                       
    //Freed memory
    subPixelStarListOnCCDResized.free();
    starListOnCCDResized.free();
    

}
//==============================================================================




//==============================================================================
/**
 * Calculates the coordinates of star and gnomonic projection into the ccd plane.
 * Compute the projection to the rotated focal plane oriented to the direction raOpticalAxis, decOpticalAxis
 * in the transformed coordinates, the origin is at raOpticalAxis, decOpticalAxis
 * @params raS  R.A. of star in degrees
 * @params decS Declination of star in degrees
 **/
void StepJitter::StepJitterTransformationSkyToCCD(double raS, double decS, double rotationAngleOA, double &xStarCCD, double &yStarCCD)
{
             

	MathTools::getGnomonicProjection(raS * Constants::DEG2RAD, decS * Constants::DEG2RAD, raOpticalAxis, decOpticalAxis, rotationAngleOA, xStarCCD, yStarCCD);
	/*	LogManager::log << "Gnomonic Projection" << std::endl;
	 LogManager::log << "RA Star raS: " << raS << std::endl;
	 LogManager::log << "Dec Star decS: " << decS << std::endl;
	 LogManager::log << "RA Projection Center raOpticalAxis: " << raOpticalAxis * Constants::RAD2DEG << std::endl;
	 LogManager::log << "Dec Projection Center decOpticalAxis: " << decOpticalAxis * Constants::RAD2DEG << std::endl;
	 LogManager::log << "Rotation Angle rotationAngleOA: " << rotationAngleOA * Constants::RAD2DEG << std::endl;
	 LogManager::log << "X Projected Star xStarCCD: " << xStarCCD << std::endl;
	 LogManager::log << "Y Projected Star yStarCCD: " << yStarCCD << std::endl;
	 */
	//the CCD is located in the focal plane with a certain offset.
	//compute the translation
	//double xCCDTrans, yCCDTrans;
    
	MathTools::translate(xStarCCD, yStarCCD, originOffsetX, originOffsetY);
	/*	LogManager::log << "Translation of CCD" << std::endl;
	 LogManager::log << "dx Translation originOffsetX: " << originOffsetX << std::endl;
	 LogManager::log << "dy Translation originOffsetY: " << originOffsetY << std::endl;
	 LogManager::log << "X Projected Translated Star xStarCCD: " << xStarCCD << std::endl;
	 LogManager::log << "Y Projected Translated Star yStarCCD: " << yStarCCD << std::endl;
	 */
	//then rotate the CCD. Rotate NOT around the origin and also not around the offset coordinates BUT around the coordinates of the original offset
	//Remember: the offset is modified in simulator.cpp because the size of the CCD is increased to correctly compute the flux at the edges
	double xRot, yRot;
	double shift = edgePixels * pixelScale * Constants::DEG2RAD / 3600.;
	MathTools::rotate(xStarCCD, yStarCCD, shift, shift, ccdOrientation * Constants::DEG2RAD, xRot, yRot);
	/*	LogManager::log << "Rotation of CCD" << std::endl;
	 LogManager::log << "Rotation center (shift): " << shift << std::endl;
	 LogManager::log << "CCD Rotation Angle ccdOrientation: " << ccdOrientation << std::endl;
	 LogManager::log << "X Projected Translated Rotated Star xRot: " << xRot << std::endl;
	 LogManager::log << "Y Projected Translated Rotated Star yRot: " << yRot << std::endl;
	 */
	double pixScalexDeg = Constants::RAD2DEG * 3600. / pixelScale;
	xStarCCD = xRot * pixScalexDeg;
	yStarCCD = yRot * pixScalexDeg;
	/*	LogManager::log << "Conversion of dimensionless (x,y)-coordinates to CCD coordinates considering the pixel scale." << std::endl;
	 LogManager::log << "Pixel scale (pixel/rad) pixScalexDeg: " << pixScalexDeg << std::endl;
	 LogManager::log << "X Pixel Coordinates Star xStarCCD: " << xStarCCD << std::endl;
	 LogManager::log << "Y Pixel Coordinate Star yStarCCD: " << yStarCCD << std::endl;
	 */
}
//==============================================================================



//==============================================================================
/**
 * Fill in the the sub-pixel map with the calculated flux values taken from the input magnitude of each detected star
 * stored in the subPixelStarListOnCCD for a given time.
 * @param myExposure Exposure time (in seconds)
 */
void StepJitter::StepJitterfillsubPixelMap(double myExposure)
{
       
    //Calculating the flux in each subpixel
    long double flux=0.0;
    double fluxMultiplicationFactor;    
    fluxMultiplicationFactor = myExposure * fluxm0 * areaTelescope * transEff * quantEff;    
    
    
    //Set the flux of each star, in its X-Y position on the subpixelMap
    for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
    {        
        //Calculate the flux detected for each star from its input magnitude and the exposure time 
        flux = fluxMultiplicationFactor * pow(10.0,  -0.4 * subPixelStarListOnCCD(i, 2));

        subPixelMap(int(round(subPixelStarListOnCCD(i, 0))), int(round(subPixelStarListOnCCD(i, 1)))) += flux;

    }
      
}
//==============================================================================


