///////////////////////////////////////////////////////////
//  DataSet.cpp
//  Implementation of the Class DataSet
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





#include "DataSet.h"




//==============================================================================
/**
* Constructor method
*/
DataSet::DataSet(){}
//==============================================================================
//==============================================================================
/**
* Destructor method
*/
DataSet::~DataSet()
{
    //pixelMap.free();             
//        subPixelMap.free();          
//        smearingMap.free();          
//        psfMap.free();               
//        cteMap.free();               
//        flatfieldMap.free();         
//        subFlatFieldMap.free();      
//        biasRegisterMap.free();      
//        jitterInputParams.free();    
//        starListOnCCD.free();        
//        subPixelStarListOnCCD.free();

}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================

unsigned int DataSet::seedRNG;


//==============================================================================
/**
*  This function uses the tinyXML library to read the input parameters file.
* 
* @param parameterFile XML input parameters file name.
*/
void DataSet::datasetReadParameterFile(string parameterFile)
{

    //Using ticpp library to read the parameters file.     
    //'Document' is a tinyXML class.
    Document doc(parameterFile);
    doc.LoadFile();

    // 'Element' is a class of tinyXML.
    Element* config = doc.FirstChildElement("Config");        
    Element* pElem = config->FirstChildElement("GeneralParameters");
    pElem->FirstChildElement("Processors")->GetText(&threads);
    pElem->FirstChildElement("MemoryLimit")->GetText(&memoryLimit);
    pElem->FirstChildElement("Convolution")->GetText(&convolutionMethod);
    pElem->FirstChildElement("UseFFTWisdom")->GetText(&useFFTWisdom);
    pElem->FirstChildElement("FFTWisdomPath")->GetText(&fftWisdomPath);
    
    pElem->FirstChildElement("SeedRNG")->GetText(&seedRNG);

    pElem = config->FirstChildElement("OutputParameters");
    pElem->FirstChildElement("OutputPath")->GetText(&outputPath);

    pElem->FirstChildElement("Prefix")->GetText(&prefix);
    pElem->FirstChildElement("WriteSubPixelFits")->GetText(&writeSubPixelFits);

    pElem = config->FirstChildElement("ObservingParameters");
    pElem->FirstChildElement("ExposureTime")->GetText(&exposureTime);
    pElem->FirstChildElement("NumExposures")->GetText(&numExposures);
    pElem->FirstChildElement("CatalogueFileName")->GetText(&catalogueFileName);

    pElem->FirstChildElement("Fluxm0")->GetText(&fluxm0);
    pElem->FirstChildElement("LightCollectingArea")->GetText(&areaTelescope);
    pElem->FirstChildElement("TransmissionEfficiency")->GetText(&transEff); 
    pElem->FirstChildElement("OpticalAxisRACenter")->GetText(&opticalAxisRACenter);
    pElem->FirstChildElement("OpticalAxisDecCenter")->GetText(&opticalAxisDecCenter);
    pElem->FirstChildElement("FocalPlaneOrientation")->GetText(&orientationFocalPlane);
    pElem->FirstChildElement("CCDPredefinedPosition")->GetText(&ccdPredefinedPosition);
    pElem->FirstChildElement("CCDOriginOffsetX")->GetText(&originOffsetXmm);
    pElem->FirstChildElement("CCDOriginOffsetY")->GetText(&originOffsetYmm);
    pElem->FirstChildElement("CCDOrientation")->GetText(&ccdOrientation);

    pElem = config->FirstChildElement("CCD");
    pElem->FirstChildElement("CCDSizeX")->GetText(&ccdSizeX);
    pElem->FirstChildElement("CCDSizeY")->GetText(&ccdSizeY);
    pElem->FirstChildElement("PixelSize")->GetText(&pixelSize);
    pElem->FirstChildElement("PixelScale")->GetText(&pixelScale);
    pElem->FirstChildElement("Gain")->GetText(&gain);
    pElem->FirstChildElement("QuantumEfficiency")->GetText(&quantEff); 
    pElem->FirstChildElement("FullWellSaturation")->GetText(&fullWellSat);
    pElem->FirstChildElement("DigitalSaturation")->GetText(&digitalSat);
    pElem->FirstChildElement("ReadoutNoise")->GetText(&readOutNoise);
    pElem->FirstChildElement("ElectronicOffset")->GetText(&electronicOffset);
    pElem->FirstChildElement("BiasPrescanRows")->GetText(&numPrescanRows);
    pElem->FirstChildElement("ReadOutTime")->GetText(&readOutTime);
    pElem->FirstChildElement("SmearingOverscanRows")->GetText(&numSmearingOverscanRows);
    pElem->FirstChildElement("FlatfieldPtPNoise")->GetText(&flatfieldPixelNoise);
    pElem->FirstChildElement("FlatfieldSubpixelNoise")->GetText(&flatfieldWhiteNoise);
    pElem->FirstChildElement("FlatfieldIntraWidth")->GetText(&flatfieldIntraPixelWidth);
    pElem->FirstChildElement("CTEMean")->GetText(&meanCTE);
    pElem->FirstChildElement("CTELowPixels")->GetText(&numLowCTEPixels);
    pElem->FirstChildElement("CTELowLines")->GetText(&numLowCTELines);

    pElem = config->FirstChildElement("SubField");
    pElem->FirstChildElement("SubFieldZeroPointX")->GetText(&subFieldZeroX);
    pElem->FirstChildElement("SubFieldZeroPointY")->GetText(&subFieldZeroY);
    pElem->FirstChildElement("SubFieldSizeX")->GetText(&subFieldSizeX);
    pElem->FirstChildElement("SubFieldSizeY")->GetText(&subFieldSizeY);
    pElem->FirstChildElement("SubPixels")->GetText(&subPixelsPerPixel);

    pElem = config->FirstChildElement("JitterParameters");
    pElem->FirstChildElement("UseJitter")->GetText(&useJitter);
    if (useJitter)
    {
        pElem->FirstChildElement("UseJitterFromFile")->GetText(&useJitterFromFile);
        pElem->FirstChildElement("JitterInterval")->GetText(&jitterInterval);
        pElem->FirstChildElement("JitterRms")->GetText(&jitterRms);
        pElem->FirstChildElement("JitterDrift")->GetText(&jitterDrift);
        pElem->FirstChildElement("JitterRepointing")->GetText(&jitterRepointing);

        pElem->FirstChildElement("JitterFileName")->GetText(&jitterFile);
        pElem->FirstChildElement("JitterAxisAngularDistance")->GetText(&jitterAngularDist);
        pElem->FirstChildElement("JitterAxisPositionAngle")->GetText(&jitterPosAngle);
        pElem->FirstChildElement("JitterMultiplicationFactor")->GetText(&jitterMultFactor);
    }

    pElem = config->FirstChildElement("PSFParameters");
    
    pElem->FirstChildElement("UseGauss")->GetText(&useGauss);
    pElem->FirstChildElement("PSFLocationDependent")->GetText(&psfLocationDependent);
    if (useGauss)
    {        
        pElem->FirstChildElement("PSFGaussFWHM")->GetText(&psfGaussFWHM);
    } 
    else if (psfLocationDependent)
    {        
        pElem->FirstChildElement("PSFLocationFileName")->GetText(&psfLocationFile);
    }         
    else
    {        
        pElem->FirstChildElement("PSFFileName")->GetText(&psfFileName);
    }         


    pElem->FirstChildElement("PSFNumRows")->GetText(&psfNumPixels);
    pElem->FirstChildElement("PSFSubPixels")->GetText(&psfSubPixels);
    pElem->FirstChildElement("PSFCenterX")->GetText(&psfCenterX);
    pElem->FirstChildElement("PSFCenterY")->GetText(&psfCenterY);
    pElem->FirstChildElement("PSFOrientation")->GetText(&psfOrientation);

    pElem->FirstChildElement("PSFRotation")->GetText(&psfRotation);
    if (psfRotation)
    {        
        pElem->FirstChildElement("PSFRotationAngle")->GetText(&psfRotationAngle);
    } 

    pElem = config->FirstChildElement("NoiseParameters");
    pElem->FirstChildElement("PhotonNoise")->GetText(&usePhotonNoise);
    pElem->FirstChildElement("SkyBackground")->GetText(&background);
    pElem->FirstChildElement("CosmicHitRate")->GetText(&cosmicHitRate);
    pElem->FirstChildElement("CosmicSaturation")->GetText(&cosmicsSatFactor);
    pElem->FirstChildElement("CosmicsLength")->GetText(&cosmicsLength);
    pElem->FirstChildElement("CosmicsWidth")->GetText(&cosmicsWidth);

    LogManager::log << "    Successfully read parameter file";
    GlobalVariables::logManager.LogManagerShowLog();        



}
//==============================================================================



//==============================================================================
/**
* This function sets the predefined CCD parameters.
* 
* @param originOffsetXmm
* @param originOffsetYmm
* @param ccdOrientation
* @param ccdSizeX
* @param ccdSizeY
* @param pixelSize
*/
void DataSet::datasetSetPredefinedCCDParams(int originOffsetXmm, int originOffsetYmm, 
    int ccdOrientation, int ccdSizeX, int ccdSizeY, int pixelSize)
{
    this->originOffsetXmm = originOffsetXmm;
    this->originOffsetYmm = originOffsetYmm;
    this->ccdOrientation = ccdOrientation;
    this->ccdSizeX = ccdSizeX;
    this->ccdSizeY = ccdSizeY;
    this->pixelSize = pixelSize;

}
//==============================================================================



//==============================================================================
/**
* This function sets the CCD parameters.
* 
* @param ccdSizeX
* @param ccdSizeY
* @param subFieldSizeX
* @param subFieldSizeY
* @param originOffsetXmm
* @param originOffsetYmm
* @param raCenterSubField
* @param declCenterSubField
* @param radiusFOVCCD
* @param edgePixels
* @param raOpticalAxis
* @param decOpticalAxis
* @param rotationAngleOA
* @param xFOVSubField
* @param yFOVSubField
*/
void DataSet::datasetSetCCDParams(int ccdSizeX, int ccdSizeY, int subFieldSizeX, int subFieldSizeY, 
            double originOffsetXmm, double originOffsetYmm,double raCenterSubField, double declCenterSubField, 
            double radiusFOVCCD, int edgePixels, double raOpticalAxis, double decOpticalAxis, double rotationAngleOA,
            double xFOVSubField, double yFOVSubField)
{
    
    this->originOffsetXmm = originOffsetXmm;
    this->originOffsetYmm = originOffsetYmm;
    this->ccdOrientation = ccdOrientation;        
    this->ccdSizeX = ccdSizeX;
    this->ccdSizeY = ccdSizeY;
    this->pixelSize = pixelSize;
    this->subFieldSizeX = subFieldSizeX;
    this->subFieldSizeY = subFieldSizeY;
    this->raCenterSubField = raCenterSubField;
    this->declCenterSubField = declCenterSubField;
    this->radiusFOVCCD = radiusFOVCCD;
    this->edgePixels = edgePixels;
    this->raOpticalAxis = raOpticalAxis; 
    this->decOpticalAxis = decOpticalAxis;         
    this->rotationAngleOA = rotationAngleOA; 
    this->xFOVSubField = xFOVSubField;
    this->yFOVSubField = yFOVSubField;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the new calculated background parameter.
* 
* @param newBackground Background flux in e-/pixel/s
*/
void DataSet::datasetSetbackground(float newBackground)
{
    this->background = newBackground;

}
//==============================================================================



//==============================================================================
/**
* x,y position of the center of the focal plane. The PSF will point towards this point.
* 
* @param xOpticalAxis x position of the center of the focal plane. 
* @param yOpticalAxis y position of the center of the focal plane.
*/
void DataSet::datasetSetCenterFocalPlane(double xOpticalAxis, double yOpticalAxis)
{
    this->xOpticalAxis = xOpticalAxis;
    this->yOpticalAxis = yOpticalAxis;

}
//==============================================================================



//==============================================================================
/**
* This function sets the PSF map in the psfMap Blitz array.
* 
* @param psfMapCopy psfMapCopy is a Blitz 2-D array with the PSF map. 
*/
void DataSet::datasetSetPSFMap(Array<float, 2> psfMapCopy)
{        
    psfMap.reference(psfMapCopy);
    this->psfMap = psfMap;
}   
//==============================================================================
//==============================================================================
///**
//* This function retrieves the PSF map.
//* 
//* @param retrievedPSFMap retrievedPSFMap addresses to the Blitz 2-D array with the PSF map. 
//*/
//void DataSet::datasetGetPSFMap(Array<float, 2> &retrievedPSFMap) 
//{        
//    //Copying the psfMap Array to retrievedPSFMap
//    retrievedPSFMap.reference(psfMap);
//}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the PSF map.
* 
*/
Array<float, 2> DataSet::datasetGetPSFMap() 
{        
    //Returning the psfMap Array 
    return psfMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the CTE map in the cteMap Blitz array.
* 
* @param cteMapCopy cteMapCopy is a Blitz 2-D array with the CTE map. 
*/
void DataSet::datasetSetCTEMap(Array<double, 2> cteMapCopy)
{        
    cteMap.reference(cteMapCopy);
    this->cteMap = cteMap;
}   
//==============================================================================
////==============================================================================
///**
//* This function retrieves the CTE map.
//* 
//* @param retrievedCTEMap retrievedCTEMap addresses to the Blitz 2-D array with the CTE map. 
//*/
//void DataSet::datasetGetCTEMap(Array<double, 2> &retrievedCTEMap) 
//{        
//    //Copying the cteMap Array to retrievedCTEMap
//    retrievedCTEMap.reference(cteMap);
//}   
////==============================================================================
//==============================================================================
/**
* This function retrieves the CTE map.
* 
*/
Array<double, 2> DataSet::datasetGetCTEMap() 
{        
    //Returning the CTE map Array 
    return cteMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the jitter displacements in the jitterInputParams Blitz array.
* 
* @param jitterInputParamsCopy jitterInputParamsCopy is a Blitz 2-D array with the jitter displacement 
*        time (in seconds), yaw, pitch and roll angles (in degrees) of the jitter time-series from the jitter input file. 
*/
void DataSet::datasetSetJitterInputParams(Array<double, 2> jitterInputParamsCopy)
{        
    jitterInputParams.reference(jitterInputParamsCopy);
    this->jitterInputParams = jitterInputParams;
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the jitter Input Params.
* 
*/
Array<double, 2> DataSet::datasetGetJitterInputParams() 
{        
    //Returning the CTE map Array 
    return jitterInputParams;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the FlatField map in the flatfieldMap Blitz array.
* 
* @param flatfieldMapCopy flatfieldMapCopy is a Blitz 2-D array with the FlatField map. 
*/
void DataSet::datasetSetFlatFieldMap(Array<float, 2> flatfieldMapCopy)
{        
    flatfieldMap.reference(flatfieldMapCopy);
    this->flatfieldMap = flatfieldMap;
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the FlatField map.
* 
* @param retrievedFlatFieldMap retrievedFlatFieldMap addresses to the Blitz 2-D array with the FlatField map. 
*/
//void DataSet::datasetGetFlatFieldMap(Array<float, 2> &retrievedFlatFieldMap) 
//{        
//    //Copying the flatfieldMap Array to retrievedFlatFieldMap
//    retrievedFlatFieldMap.reference(flatfieldMap);
//}   
//==============================================================================


//==============================================================================
/**
* This function retrieves the FlatField map.
* 
*/
Array<float, 2> DataSet::datasetGetFlatFieldMap() 
{        
    //Returning the flatfieldMap Array 
    return flatfieldMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the subFlatField map in the subFlatFieldMap Blitz array.
* 
* @param subflatfieldMapCopy subflatfieldMapCopy is a Blitz 2-D array with the subFlatField map. 
*/
void DataSet::datasetSetsubFlatFieldMap(Array<float, 2> subflatfieldMapCopy)
{        
    subFlatFieldMap.reference(subflatfieldMapCopy);
    this->subFlatFieldMap = subFlatFieldMap;
}   
//==============================================================================
////==============================================================================
///**
//* This function retrieves the subFlatField map.
//* 
//* @param retrievedsubFlatFieldMap retrievedsubFlatFieldMap addresses to the Blitz 2-D array with the subFlatField map. 
//*/
//void DataSet::datasetGetsubFlatFieldMap(Array<float, 2> &retrievedsubFlatFieldMap) 
//{        
//    //Copying the subFlatFieldMap Array to retrievedsubFlatFieldMap
//    retrievedsubFlatFieldMap.reference(subFlatFieldMap);
//}   
////==============================================================================
//==============================================================================
/**
* This function retrieves the subFlatField map.
* 
*/
Array<float, 2> DataSet::datasetGetsubFlatFieldMap() 
{        
    //Returning the subFlatFieldMap Array 
    return subFlatFieldMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the pixel maps.
* 
* @param pixelMapCopy pixelMapCopy is a Blitz 2-D array with the pixel map. 
*/
void DataSet::datasetSetpixelMap(Array<float, 2> pixelMapCopy)
{
    pixelMap.reference(pixelMapCopy);
    this->pixelMap = pixelMap;        
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the pixel map.
* 
*/
Array<float, 2>  DataSet::datasetGetpixelMap() 
{        
    //Returning the PixelMap Array 
    return pixelMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the subpixel map.
* 
* @param subPixelMapCopy subPixelMapCopy is a Blitz 2-D array with the subpixel map. 
*/
void DataSet::datasetSetsubPixelMap(Array<float, 2> subPixelMapCopy)
{
    subPixelMap.reference(subPixelMapCopy);
    this->subPixelMap = subPixelMap;        
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the subpixel map.
* 
*/
Array<float, 2>  DataSet::datasetGetsubPixelMap() 
{        
    //Returning the subPixelMap Array 
    return subPixelMap;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the initial pixel maps.
* 
* @param initPixelMapCopy initPixelMapCopy is a Blitz 2-D array with the pixel map. 
*/
void DataSet::datasetSetinitPixelMap(Array<float, 2> initPixelMapCopy)
{
    initPixelMap.reference(initPixelMapCopy);
    this->initPixelMap = initPixelMap;      
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the initial pixel map.
* 
*/
Array<float, 2> DataSet::datasetGetinitPixelMap() 
{        
    //Returning the initPixelMap 
    return  initPixelMap;
} 
//==============================================================================



//==============================================================================
/**
* This function sets the initial subpixel map.
* 
* @param initSubPixelMapCopy initSubPixelMapCopy is a Blitz 2-D array with the subpixel map. 
*/
void DataSet::datasetSetinitSubPixelMap(Array<float, 2> initSubPixelMapCopy)
{
    initSubPixelMap.reference(initSubPixelMapCopy);
    this->initSubPixelMap = initSubPixelMap;        
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the initial subpixel map.
* 
*/
Array<float, 2> DataSet::datasetGetinitSubPixelMap() 
{        
    //Returning the initSubPixelMap 
    return  initSubPixelMap;
} 
//==============================================================================



//==============================================================================
/**
* This function sets the smearing map.
* 
* @param smearingMapCopy smearingMapCopy is a Blitz 2-D array with the shielded register 
*        of the CCD that only contains the star smearing due to frame transfer. 
*/
void DataSet::datasetSetsmearingMap(Array<float, 2> smearingMapCopy)
{
    smearingMap.reference(smearingMapCopy);
    this->smearingMap = smearingMap;        
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the smearing map at a Blitz 2-D array with the shielded 
 * register of the CCD that only contains the star smearing due to frame transfer. 
*        
*/
Array<float, 2> DataSet::datasetGetsmearingMap() 
{        
    //Returning the smearingMap Array 
    return smearingMap;
}   
//==============================================================================




//==============================================================================
/**
* This function sets the bias overscan strip of the CCD map.
* 
* @param biasRegisterMapCopy biasRegisterMapCopy is a Blitz 2-D array with the bias 
* overscan strip of the CCD map. 
*/
void DataSet::datasetSetbiasRegisterMap(Array<float, 2> biasRegisterMapCopy)
{
    biasRegisterMap.reference(biasRegisterMapCopy);
    this->biasRegisterMap = biasRegisterMap;        
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the bias overscan strip of the CCD map.
*        
*/
Array<float, 2> DataSet::datasetGetbiasRegisterMap() 
{        
    //Returning the biasRegisterMap Array 
    return biasRegisterMap;
}   
//==============================================================================




//==============================================================================
/**
* This function sets the list of stars on the CCD in the starListOnCCD Blitz array.
* 
* @param starListOnCCDCopy starListOnCCDCopy is a Blitz 2-D array with the list of stars on the CCD. 
*/
void DataSet::datasetSetStarListOnCCD(Array<float, 2> starListOnCCDCopy)
{        

    starListOnCCD.reference(starListOnCCDCopy);
    this->starListOnCCD = starListOnCCD;
}   
//==============================================================================
//==============================================================================
/**
* This function retrieves the list of stars on the CCD in the starListOnCCD map.
*        
*/
Array<float, 2> DataSet::datasetGetStarListOnCCD() 
{        
    //Returning the starListOnCCD Array 
    return starListOnCCD;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the list of stars on the CCD in the subPixelStarListOnCCD Blitz array.
* 
* @param subPixelStarListOnCCDCopy subPixelStarListOnCCDCopy is a Blitz 2-D array with the list of stars on the CCD. 
*/
void DataSet::datasetSetsubPixelStarListOnCCD(Array<float, 2> subPixelStarListOnCCDCopy)
{        

    subPixelStarListOnCCD.reference(subPixelStarListOnCCDCopy);
    this->subPixelStarListOnCCD = subPixelStarListOnCCD;

}
//==============================================================================
//==============================================================================
/**
* This function retrieves the list of stars on the CCD in the subPixelStarListOnCCD map.
* 
*/      
Array<float, 2> DataSet::datasetGetsubPixelStarListOnCCD() 
{        
    //Returning the subPixelStarListOnCCD Array 
    return subPixelStarListOnCCD;
}   
//==============================================================================



//==============================================================================
/**
* This function sets the array with the exposures names in the exposuresNamesArray Blitz array.
* This names will be used for naming the exposures output files.
*
* @param exposuresNamesArrayCopy exposuresNamesArrayCopy is a Blitz 2-D array with the list of exposures names.
*/
void DataSet::datasetSetExposuresNamesArray(Array<string, 1> exposuresNamesArrayCopy)
{

    exposuresNamesArray.reference(exposuresNamesArrayCopy);
    this->exposuresNamesArray = exposuresNamesArray;
}
//==============================================================================
//==============================================================================
/**
* This function retrieves the list of stars of exposures names in the exposuresNamesArray.
*
*/
Array<string, 1> DataSet::datasetGetExposuresNamesArray()
{
    //Return the exposuresNamesArray Array 
    return exposuresNamesArray;
}
//==============================================================================





//==============================================================================
/**
* This function sets the star catalogue in the starCatalogue Blitz array.
* 
* @param starCatalogueCopy starCatalogueCopy is a Blitz 2-D array with the star catalogue
*        containing id, RA, declination and magnitude of each star in the input star catalogue.
*/
void DataSet::datasetSetStarCatalogue(Array<float, 2> starCatalogueCopy)
{

    starCatalogue.reference(starCatalogueCopy);
    this->starCatalogue = starCatalogue;
}
//==============================================================================
//==============================================================================
/**
* This function retrieves the list of stars  in the starCatalogue map in a Blitz 2-D array 
* with the star catalogue containing id, RA, declination and magnitude of each star 
* in the input star catalogue.
*/
Array<float, 2> DataSet::datasetGetStarCatalogue()
{
    //Returning the starCatalogue Array 
    return starCatalogue;
}
//==============================================================================






//==============================================================================
/**
* Constructor method
*/
DataSetPhotometry::DataSetPhotometry(){}
//==============================================================================
//==============================================================================
/**
* Destructor method
*/
DataSetPhotometry::~DataSetPhotometry(){}
//==============================================================================




//==============================================================================
/**
*  This function uses the tinyXML library to read the photometry input parameters file.
* 
* @param photometryParameterFile Photometry input parameters file name.
*/
void DataSetPhotometry::datasetPhotometryReadParameterFile(string photometryParameterFile)
{

    //Using ticpp library to read the photometry parameters file.
    LogManager::log << "    Reading photometry parameter file ";
    GlobalVariables::logManager.LogManagerShowLog();


    //'Document' is a tinyXML class.
    Document photometryFile(photometryParameterFile);
    photometryFile.LoadFile();
    LogManager::log << "    Successfully opened parameter file";
    GlobalVariables::logManager.LogManagerShowLog();

    Element* pElem = photometryFile.FirstChildElement ( "GeneralParameters" );
    pElem->FirstChildElement ( "PhotometryMethod" )->GetText ( &photometryMethod );
    pElem->FirstChildElement ( "NumberTelescopes" )->GetText ( &photometryNumTelescopes );
    pElem->FirstChildElement ( "CorrectFlatfield" )->GetText ( &flatfieldCorrection );
    pElem->FirstChildElement ( "CorrectTrails" )->GetText ( &frameTransferSmearingCorrection );

    pElem = photometryFile.FirstChildElement ( "BackgroundParameters" );
    pElem->FirstChildElement ( "BackgroundFlux" )->GetText ( &photometryBackground );
    pElem->FirstChildElement ( "BackgroundAnnulusInnerRadius" )->GetText ( &backgroundAnnulusInnerRadius );
    pElem->FirstChildElement ( "BackgroundAnnulusOuterRadius" )->GetText ( &backgroundAnnulusOuterRadius );


    if ( photometryMethod == "WM" ) //if Weighted Mask.
    {
            photometryDirName += "photometry_weighted";
    }
    else if ( photometryMethod == "AP" )    //if Simple Aperture Photometry
    {
            photometryDirName += "photometry_aperture";
    }


}
//==============================================================================




//==============================================================================
/**
* This function sets the photometry Parameters.
* 
* @param photometryDirName photometryDirName is output directory for photometry. 
* @param photometryPlotsDir photometryPlotsDir is is output directory for the photometry plots. 
*/
void DataSetPhotometry::datasetSetPhotometryParams(string photometryDirName, string photometryPlotsDir)
{        

    this->photometryDirName = photometryDirName;
    this->photometryPlotsDir = photometryPlotsDir;

}   
//==============================================================================



