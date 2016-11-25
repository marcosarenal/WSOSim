///////////////////////////////////////////////////////////
//  ParamsStarfield.cpp
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



#include "ParamsStarfield.h"


//==============================================================================
/**
* Constructor method
*/
ParamsStarfield::ParamsStarfield(){}
//==============================================================================
//==============================================================================
/**
* Destructor method
*/
ParamsStarfield::~ParamsStarfield(){}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================
/**
* Define a new starfield object which contains the x,y-coordinates on the detector
* of the stars in the detector frame.
*/
void ParamsStarfield::ParamsStarfieldSetting(DataSet &m_DataSet)

{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    //Getting Catalogue File Name from DataSet
    catalogueFileName = p_DataSet->datasetGetCatalogueFileName();

    //Getting detector parameters
    pixelScale = p_DataSet->datasetGetpixelScale();
    originOffsetXmm = p_DataSet->datasetGetOriginOffsetXmm();
    originOffsetYmm = p_DataSet->datasetGetOriginOffsetYmm();

    pixelSize = p_DataSet->datasetGetpixelSize();
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();

    originOffsetX = originOffsetXmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);
    originOffsetY = originOffsetYmm * pixelScale * Constants::DEG2RAD / (pixelSize * 3.6);

    edgePixels = p_DataSet->datasetGetedgePixels();

    //Call the ParamsStarfieldReadCatalogue function to read the Starfield Catalogue
    //and find stars within the circular FOV 
    ParamsStarfield::ParamsStarfieldReadCatalogue(catalogueFileName);

    //Compute the x,y-coordinates on the detector of the stars in the starfield
    ParamsStarfield::ParamsStarfieldCoords(m_DataSet);

    //Fill in the pixel map and the subpixel map
    ParamsStarfield::ParamsStarfieldfillsubPixelMap(m_DataSet);
    
   
    //Compute the x,y position of the center of the focal plane. The PSF will point towards this point
    ParamsStarfield::ParamsStarfieldTransformationSkyToCCD(Constants::RAD2DEG * raOpticalAxis, Constants::RAD2DEG * decOpticalAxis, xOpticalAxis, yOpticalAxis);

    //Set this xOpticalAxis and yOpticalAxis parameters in the DataSet to be used in the StepConvolvePSF module.
    p_DataSet->datasetSetCenterFocalPlane(xOpticalAxis, yOpticalAxis);




    LogManager::log << "    Found " << numStarsCCD << " stars on CCD and " << numStarsSubField
    << " stars on sub-image.";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

    //TODO: update these outputs (pma).
    //	LogManager::log << "Center of sub-image field is at: RA=" << ccd.getRACenterSubField() << ", DEC="
    //			<< ccd.getDeclCenterSubField();
    //	GlobalVariables::logManager.LogManagerAppendLogAndShow();
    //	LogManager::log << "Field of view of sub-image in x-direction is: " << ccd.getXFOVSubField();
    //	GlobalVariables::logManager.LogManagerAppendLogAndShow();
    //	LogManager::log << "Field of view of sub-image in y-direction is: " << ccd.getYFOVSubField();
    //	GlobalVariables::logManager.LogManagerAppendLogAndShow();
    //	LogManager::log << "Orientation of the jitter roll axis: RA=" << raJitterAxis << " , DEC=" << declJitterAxis;
    //	double starsperdeg = ccd.getNumStarsSubField() / (ccd.getXFOVSubField() * ccd.getYFOVSubField() * Constants::Pi);
    //	GlobalVariables::logManager.LogManagerAppendLogAndShow();
    //	LogManager::log << "Approximate number of stars per square degree on sub-field: " << starsperdeg;
    //	GlobalVariables::logManager.LogManagerAppendLogAndShow();

    //Check whether there are stars in the field to continue the simulation or exit.
    if (numStarsSubField == 0)
    {
        std::string input = "";
        LogManager::log << "\nWARNING: No stars on sub-field! Continue? (y/n)" << std::endl;
        GlobalVariables::logManager.LogManagerShowLog();
        getline(cin, input);
        if (input != "y" && input != "Y")
        {
            LogManager::log << "No stars on sub-field. Simulation terminated. ";
            GlobalVariables::logManager.LogManagerAppendLogAndShow();
            exit(1);
        }
    }

}
//==============================================================================




//==============================================================================
/**
* Read a star catalog and write all stars to an array that are inside the given
* radius centered on (centerRA, centerDec). Since a circular area is read in,
* the arrays will contain more stars than are actually on the CCD.
* @param catalogueFileName File name of star catalog. The catalog must be an ASCII
* file consisting of 3 columns separated by spaces: R.A, Dec and magnitude.
*/
void ParamsStarfield::ParamsStarfieldReadCatalogue(std::string catalogueFileName)
{
    //read in the star catalouge
    LogManager::log<< "    Reading the star catalogue from " << catalogueFileName;
    GlobalVariables::logManager.LogManagerAppendLogAndShow();


    //Check whether the catalogue file can be opened or not
    ifstream myfile(catalogueFileName.c_str());
    if (!myfile.is_open())
    {
        LogManager::log<< "Error: Unable to open file " << catalogueFileName;
        GlobalVariables::logManager.LogManagerAppendLogAndShow();

        exit(0);
    }

    //There must be one line in the catalogue file for each star
    int numLines = FileUtilities::countLines(catalogueFileName);

    if (numLines == 0)
    {
        LogManager::log<< "Error: Empty file " << catalogueFileName << ".";
        GlobalVariables::logManager.LogManagerAppendLogAndShow();

        exit(0);
    }


    //Retrieving parameters from DataSet
    //Radius of the FOV in degrees (should at least cover the complete CCD).
    radiusFOVCCD = p_DataSet->datasetGetradiusFOVCCD();

    //Right Ascension of center of FOV of CCD in degrees.
    raCenterSubField = p_DataSet->datasetGetraCenterSubField();

    //Declination of center of FOV of CCD in degrees.
    declCenterSubField = p_DataSet->datasetGetdeclCenterSubField();


    //Initialization of the ra, dec, id and magn parameters
    ra.resize(numLines);
    ra = 0.0;
    
    dec.resize(numLines);
    dec = 0.0;
    
    magn.resize(numLines);
    magn = 0.0;
    
    id.resize(numLines);
    id = 0;

    int numb_stars = 0;
    int i = 0;

    double myRA;
    double myDec;
    double myMag;

    //From degrees to radians.
    double raCenterRad = raCenterSubField * Constants::DEG2RAD;
    double decCenterRad = declCenterSubField * Constants::DEG2RAD;
    double radiusFOV = (1. + radiusFOVCCD) * Constants::DEG2RAD;

    if (radiusFOV > Constants::Pid2)
    {
        radiusFOV = Constants::Pid2;
    }

    LogManager::log << "    Radius of Field of view of CCD is: " << radiusFOV * Constants::RAD2DEG  <<" degrees" ;
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

    LogManager::log << "    Center of CCD field is at: RA=" << raCenterRad * Constants::RAD2DEG << " degrees,"
    " DEC=" << decCenterRad* Constants::RAD2DEG<<" degrees.";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

    //Reading each line of the catalogue.
    try
    {        
        while (myfile >> myRA >> myDec >> myMag)
        {
            //if the star is within the FOV
            if (MathTools::getAngularDistanceSphere(raCenterRad, decCenterRad, myRA * Constants::DEG2RAD, myDec * Constants::DEG2RAD)
                < radiusFOV) 
            {
                ra(numb_stars) = myRA;
                dec(numb_stars) = myDec;
                magn(numb_stars) = myMag;
                id(numb_stars) = i;

                numb_stars++;
            }
            i++;
        }
    }
    catch (int e)
    {
        std::cerr << "An exception occurred. Exception Nr. " << e  << "There seems to be a problem with the star catalogue."<<std::endl;
        std::cerr << "It might not be an appropriate format. Remember: RA, DEC, magnitude (in columns).   "<<std::endl;
    }
        
    myfile.close();
    
    
    //Check whether there are stars in the field to continue the simulation or exit.
    if (numb_stars == 0)
    {
        LogManager::log << "No stars on sub-field. Simulation terminated. ";
        GlobalVariables::logManager.LogManagerAppendLogAndShow();
        exit(1);
        
    }
    

    LogManager::log<< "    Found " << numb_stars << " stars in the star catalogue. ";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

}
//==============================================================================





//==============================================================================
/**
* Setting the star list of the starfield on the CCD into the starListOnCCD array.
* The starListOnCCD contains for each star, its X and Y position in CCD, magnitude, RA,
* declination and identification number.
*/
void ParamsStarfield::ParamsStarfieldCoords(DataSet &m_DataSet)
{
    //Initialization of parameters
    NumStars = magn.rows();
    numStarsCCD = 0;
    numStarsSubField = 0;

    //Retrieving parameters
    ccdSizeX = p_DataSet->datasetGetccdSizeX();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();
    subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY = p_DataSet->datasetGetsubFieldZeroY();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    ccdOrientation = p_DataSet->datasetGetccdOrientation();
    raOpticalAxis = p_DataSet->datasetGetraOpticalAxis();
    decOpticalAxis = p_DataSet->datasetGetdecOpticalAxis();
    rotationAngleOA = p_DataSet->datasetGetrotationAngleOA();


    // Defining the starListOnCCD blitz array
    starListOnCCD.resize(NumStars, 6);
    subPixelStarListOnCCD.resize(NumStars, 4);
    
    starCatalogue.resize(NumStars, 6);
        
    double xTrue = 0, yTrue = 0;
    int xTrueInt, yTrueInt;

    //Iterating the star list to check whether it is on the CCD and calculate its position if it is.
    for (int i = 0; i < NumStars; i++)
    {        

        //Compute the x,y-coordinates on the detector of the stars in the starfield
        ParamsStarfield::ParamsStarfieldTransformationSkyToCCD(ra(i), dec(i), xTrue, yTrue);

        //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at
        //intrapixel level (because we round and int)
        xTrueInt = int(round(xTrue * subPixelsPerPixel));
        yTrueInt = int(round(yTrue * subPixelsPerPixel));

        if (xTrueInt >= 0 && yTrueInt >= 0 && xTrueInt < ccdSizeX * subPixelsPerPixel && yTrueInt < ccdSizeY * subPixelsPerPixel)
        {
            //star is in the FoV of the CCD
            starListOnCCD(numStarsCCD, 0) = xTrue * subPixelsPerPixel; //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at
            starListOnCCD(numStarsCCD, 1) = yTrue * subPixelsPerPixel; //intrapixel level (because we round and int)
            starListOnCCD(numStarsCCD, 2) = magn(i);
            starListOnCCD(numStarsCCD, 3) = ra(i);
            starListOnCCD(numStarsCCD, 4) = dec(i);
            starListOnCCD(numStarsCCD, 5) = id(i);
                
            numStarsCCD++;

            if (xTrueInt >= subFieldZeroX * subPixelsPerPixel && yTrueInt >= subFieldZeroY * subPixelsPerPixel && xTrueInt
                < (subFieldZeroX + subFieldSizeX) * subPixelsPerPixel && yTrueInt < (subFieldZeroY + subFieldSizeY)
                * subPixelsPerPixel)
            {
                //star is in the FoV of the sub image
                subPixelStarListOnCCD(numStarsSubField, 0) = (xTrue - subFieldZeroX) * subPixelsPerPixel; //multiply all with subPixelsPerPixel to ensure that we check if the star is on the frame at 
                subPixelStarListOnCCD(numStarsSubField, 1) = (yTrue - subFieldZeroY) * subPixelsPerPixel; //intrapixel level (because we round and int)
                subPixelStarListOnCCD(numStarsSubField, 2) = magn(i);
                subPixelStarListOnCCD(numStarsSubField, 3) = id(i);
            
                numStarsSubField++;
                    
            }
        }        
        
        //Copy the star catalogue to the DataSet
        starCatalogue(i, 0) = id(i);
        starCatalogue(i, 1) = ra(i);
        starCatalogue(i, 2) = dec(i);
        starCatalogue(i, 3) = magn(i);
                
    }
        

    //Set the starCatalogue into the DataSet            
    m_DataSet.datasetSetStarCatalogue(starCatalogue);

    // Defining a blitz array to remove zeros from the starListOnCCD array
    Array<float, 2> subPixelStarListOnCCDResized;
    Array<float, 2> starListOnCCDResized;

    starListOnCCDResized.resize(numStarsCCD, 6);
    starListOnCCDResized = 0.0;
    starListOnCCDResized(Range::all(),Range::all()) = starListOnCCD(Range(0,numStarsCCD), Range::all());

    starListOnCCD.resize(numStarsCCD, 6);
    starListOnCCD = 0.0;
    starListOnCCD = starListOnCCDResized;
        

    //Remove zeros from the subPixelStarListOnCCD array
    subPixelStarListOnCCDResized.resize(numStarsSubField, 4);
    subPixelStarListOnCCDResized(Range::all(), Range::all()) = subPixelStarListOnCCD(Range(0,numStarsSubField), Range::all());


    subPixelStarListOnCCD.resize(numStarsSubField, 4);
    subPixelStarListOnCCD = subPixelStarListOnCCDResized;

                    
    //Setting the star list in the DataSet to be used in StepJitter
    m_DataSet.datasetSetsubPixelStarListOnCCD(subPixelStarListOnCCD);

    //Setting the star list in the DataSet to be used in StepJitter
    m_DataSet.datasetSetStarListOnCCD(starListOnCCD);

        

    //Freed memory
    subPixelStarListOnCCDResized.free();
    starListOnCCDResized.free();
                        
    ra.free();
    dec.free();
    id.free();
    magn.free();

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
void ParamsStarfield::ParamsStarfieldTransformationSkyToCCD(double raS, double decS, double &xStarCCD, double &yStarCCD)
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
* Fill in the the sub-pixel map with the calculated flux values.
*/
void ParamsStarfield::ParamsStarfieldfillsubPixelMap(DataSet &m_DataSet)
{

    //Getting exposureTime from DataSet
    exposureTime = p_DataSet->datasetGetexposureTime();

    // Retrieving the empty sub-pixel map from DataSet
    //Initialize the initSubPixelMap
    initSubPixelMap.resize(m_DataSet.datasetGetinitSubPixelMap().extent(0),m_DataSet.datasetGetinitSubPixelMap().extent(1));
    initSubPixelMap = 0.0;
    //Retrieving the initSubPixelMap from DataSet            
    initSubPixelMap = m_DataSet.datasetGetinitSubPixelMap();

    //Retrieving the subPixelStarListOnCCD contains for each star, its X and Y position in CCD, magnitude, RA, 
    //declination and identification number.
    subPixelStarListOnCCD.resize(m_DataSet.datasetGetsubPixelStarListOnCCD().extent(0),m_DataSet.datasetGetsubPixelStarListOnCCD().extent(1));
    subPixelStarListOnCCD = 0.0;
    subPixelStarListOnCCD = m_DataSet.datasetGetsubPixelStarListOnCCD();


    // Retrieving parameters from DataSet
    fluxm0 = p_DataSet->datasetGetfluxm0();
    areaTelescope  = p_DataSet->datasetGetareaTelescope();
    transEff = p_DataSet->datasetGettransEff();
    quantEff = p_DataSet->datasetGetquantEff();

    //Calculating the flux in each subpixel
    float flux=0.0;

            
    for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
    {
        if(int(round(subPixelStarListOnCCD(i, 0)))<initSubPixelMap.rows())
        {
            flux = exposureTime * fluxm0 * pow(10.0, -0.4 * subPixelStarListOnCCD(i, 2)) * areaTelescope * transEff * quantEff;
            initSubPixelMap(int(round(subPixelStarListOnCCD(i, 0))), int(round(subPixelStarListOnCCD(i, 1)))) += flux;
        }
    }
    

    //Set the calculated initSubPixelMap in the DataSet
    m_DataSet.datasetSetinitSubPixelMap(initSubPixelMap);

}
//==============================================================================
