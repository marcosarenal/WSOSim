///////////////////////////////////////////////////////////
//  StepPhotometryOnTheGo.cpp
//  Implementation of the Class StepPhotometryOnTheGo
//  Created on:      September 24, 2013, 12:05 PM
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



#include "StepPhotometryOnTheGo.h"


//==============================================================================
/**
 * Constructor method
 */
StepPhotometryOnTheGo::StepPhotometryOnTheGo(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepPhotometryOnTheGo::~StepPhotometryOnTheGo(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function performs the photometry processing ON THE GO.
 * This Processing Step has the same functionality as the ProcessingPhotometry::StepPhotometryOnTheGoPipeline() method
 * but here is implemented differently in order to apply the photometry to each image right after it is created
 * instead of reading its FITS file image in order to save computing time.
 * The flux of each star in the 
 * sub-field is measured assuming a Gaussian weighted mask and then derived the 
 * noise-to-signal ratio (NSR).
 */
void StepPhotometryOnTheGo::StepPhotometryOnTheGoPipeline(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    p_datasetPhotometry = &datasetPhotometry;
    
    //Retrieving parameters from DataSet
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels();
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();
    gain = p_DataSet->datasetGetgain();
    numSmearingOverscanRows  = p_DataSet->datasetGetnumSmearingOverscanRows();
    numPrescanRows = p_DataSet->datasetGetnumPrescanRows();
    photometryNumTelescopes = datasetPhotometry.datasetGetphotometryNumTelescopes();
    digitalSat = p_DataSet->datasetGetdigitalSat();
   
    exposureTime = p_DataSet->datasetGetexposureTime();
    edgePixels = p_DataSet->datasetGetedgePixels();
    transEff = p_DataSet->datasetGettransEff();
    quantEff = p_DataSet->datasetGetquantEff();
    fluxm0 = p_DataSet->datasetGetfluxm0();
    areaTelescope = p_DataSet->datasetGetareaTelescope();
    readOutTime = p_DataSet->datasetGetreadOutTime();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    

    //Retrieving the exposures Names Array from DataSet
    //Initialize exposuresNamesArray      
    exposuresNamesArray.resize(m_DataSet.datasetGetExposuresNamesArray().extent(0));
    exposuresNamesArray = m_DataSet.datasetGetExposuresNamesArray();


    //Retrieving the PSF Map from DataSet
    //Initialize PSF map  
    psfMap.resize(m_DataSet.datasetGetPSFMap().extent(0),m_DataSet.datasetGetPSFMap().extent(1));
    psfMap = 0.0;
    psfMap = m_DataSet.datasetGetPSFMap();
    
    
    //Define a flux correction factor to avoid repeating operations
    fluxCorrectionFactor = fluxm0 * areaTelescope * transEff * quantEff;
    
    //Retrieving parameters from PhotometryDataSet.
    photometryDirName = datasetPhotometry.datasetGetphotometryDirName();
    photometryPlotsDir = datasetPhotometry.datasetGetphotometryPlotsDir();
    
    photometryMethod = datasetPhotometry.datasetGetphotometryMethod();
    photometryBackground = datasetPhotometry.datasetGetphotometryBackground();
    
    frameTransferSmearingCorrection = datasetPhotometry.datasetGetframeTransferSmearingCorrection();
    flatfieldCorrection = datasetPhotometry.datasetGetphotometryFlatfieldCorrectionr();
    
    backgroundAnnulusInnerRadius = datasetPhotometry.datasetGetphotometryBackgroundAnnulusInnerRadius();
    backgroundAnnulusOuterRadius = datasetPhotometry.datasetGetphotometryBackgroundAnnulusOuterRadius();
    
    
    //Initializing and retrieving the pixel map from DataSet
    pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
    pixelMap = 0.0;
    pixelMap = m_DataSet.datasetGetpixelMap();
   
    //Retrieving the subPixelStarListOnCCD contains for each star, its X and Y position in CCD, magnitude, RA, 
    //declination and identification number.
    subPixelStarListOnCCD.resize(m_DataSet.datasetGetsubPixelStarListOnCCD().extent(0),m_DataSet.datasetGetsubPixelStarListOnCCD().extent(1));
    subPixelStarListOnCCD = 0.0;
    subPixelStarListOnCCD = m_DataSet.datasetGetsubPixelStarListOnCCD();

    
    //Retrieving the BIAS register map from DataSet        
    //Initialize smearingMap
    biasRegisterMap.resize(m_DataSet.datasetGetbiasRegisterMap().extent(0), m_DataSet.datasetGetbiasRegisterMap().extent(1));
    biasRegisterMap = 0.0;
    biasRegisterMap = m_DataSet.datasetGetbiasRegisterMap();

    //Retrieving the smearing map from DataSet
    smearingMap.resize(m_DataSet.datasetGetsmearingMap().extent(0), m_DataSet.datasetGetsmearingMap().extent(1));
    smearingMap = 0.0;
    smearingMap = m_DataSet.datasetGetsmearingMap();
    
    
    
    //Compute the Weighted PSF mask.
    StepPhotometryOnTheGo::StepPhotometryOnTheGoPreComputeWeightedMask();
            
    //Set the flux values as taken from the pixelMap.
    Array<float, 2> fluxMap(pixelMap.rows(), pixelMap.cols());
                    
    fluxMap = pixelMap;
    
    //Remove trailing from flux map using the star smearing (trailing) register
    if (frameTransferSmearingCorrection && numSmearingOverscanRows > 0)
    {
        for (int j = 0; j < smearingMap.rows(); j++)
        {    
            fluxMap(j, Range::all()) -= MathTools::median(smearingMap(j, Range::all()));            
        }
    }
           
    //After subtraction of trails, flux might be negative! Set all negative values to 0
    fluxMap = where(fluxMap < 0., 0., fluxMap);


    //Correct flux with gain, exposure time TODO: and the number of telescopes (pma)
    fluxMap *= (gain / (exposureTime));
    
    //Compute the weighted mask
    StepPhotometryOnTheGo::StepPhotometryOnTheGoMeasureStarsFluxes(fluxMap);

    LogManager::log << "    Successfully performed photometry on-the-go.";
    GlobalVariables::logManager.LogManagerShowLog();
   
    
}
//==============================================================================





//==============================================================================
/**
 * This function takes the subPixelStarListOnCCD array and recompute the star coordinates on the CCD to take the cutting of the edges into account.
 */
void StepPhotometryOnTheGo::StepPhotometryOnTheGoSubPixelStarListOnCCDCut()
{
    Array<float, 2>  tempSubPixelStarListOnCCD;
	int starCounter = 0;
	int x, y;
	tempSubPixelStarListOnCCD.resize(subPixelStarListOnCCD.rows(), subPixelStarListOnCCD.cols());
    tempSubPixelStarListOnCCD = 0.0;

	for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
	{
		x = int(round(subPixelStarListOnCCD(i, 0)));
		y = int(round(subPixelStarListOnCCD(i, 1)));
		if (x >= edgePixels * subPixelsPerPixel  && x < (subFieldSizeX - edgePixels) * subPixelsPerPixel  && 
                y >= edgePixels * subPixelsPerPixel && y < (subFieldSizeY - edgePixels) * subPixelsPerPixel)
		{
			tempSubPixelStarListOnCCD(starCounter, 0) = subPixelStarListOnCCD(i, 0) - double(2*edgePixels * subPixelsPerPixel);
			tempSubPixelStarListOnCCD(starCounter, 1) = subPixelStarListOnCCD(i, 1) - double(2*edgePixels * subPixelsPerPixel);
			tempSubPixelStarListOnCCD(starCounter, 2) = subPixelStarListOnCCD(i, 2);
			tempSubPixelStarListOnCCD(starCounter, 3) = subPixelStarListOnCCD(i, 3);
            
            starCounter++;
		}
	}    

    
	tempSubPixelStarListOnCCD.resizeAndPreserve(starCounter, 4);
    
    subPixelStarListOnCCD.resize(tempSubPixelStarListOnCCD.rows(), tempSubPixelStarListOnCCD.cols());
    subPixelStarListOnCCD = 0.0;

    subPixelStarListOnCCD(Range::all(), Range::all()) = tempSubPixelStarListOnCCD(Range::all(), Range::all());
    
    
    tempSubPixelStarListOnCCD.free();
                
    //Set the resized subPixelStarListOnCCD into the DataSte to be retrieved afterwards
    p_DataSet->datasetSetsubPixelStarListOnCCD(subPixelStarListOnCCD);
    

	subFieldSizeX -= 2 * edgePixels;
	subFieldSizeY -= 2 * edgePixels;
}
//==============================================================================




//==============================================================================
/**
 * This function computes the weighted mask using the algorithm as proposed by Reza Samadi.
 */
void StepPhotometryOnTheGo::StepPhotometryOnTheGoPreComputeWeightedMask()
{
    
    //Retrieving parameters from DataSet
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels()+1;
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    photometryMethod = p_datasetPhotometry->datasetGetphotometryMethod();
    
    //Create a bigger PSF map
	Array<float, 2> biggerPSF(psfSubPixels + 2 * subPixelsPerPixel - 1, psfSubPixels + 2 * subPixelsPerPixel - 1);
    
    //Initialize the weighted mask
	weightedPSF.resize(subPixelsPerPixel + 1, subPixelsPerPixel + 1, int(biggerPSF.rows() / subPixelsPerPixel), int(biggerPSF.cols() / subPixelsPerPixel));
    weightedPSF = 0.0;
    
                
    
	normalizationFactor.resize(subPixelsPerPixel + 1, subPixelsPerPixel + 1);
    normalizationFactor = 0.0;  
   
    //For each subpixel
	for (int i = 0; i <= subPixelsPerPixel; i++)
    {
		for (int j = 0; j <= subPixelsPerPixel; j++)
		{
            
			//shift PSF to the closest subcoordinates of the star
			biggerPSF = 0.;
			biggerPSF(Range(i, i + psfMap.cols() - 1), Range(j, j + psfMap.cols() - 1)) = psfMap;
            
			//rebin biggerPSF to normal pixel scale
			for (int l = 0; l < weightedPSF.extent(thirdDim); l++)
            {
				for (int k = 0; k < weightedPSF.extent(fourthDim); k++)
                {
					weightedPSF(i, j, l, k) = sum(biggerPSF(Range(l * subPixelsPerPixel, (l + 1) * subPixelsPerPixel - 1), 
                                                            Range(k * subPixelsPerPixel, (k + 1) * subPixelsPerPixel - 1)));
                }
            }
            
            //if photometryMethod == WeightedMask
            if (photometryMethod == "WM") 
            {
                double c = sum(pow2(weightedPSF(i, j, Range::all(), Range::all())));
                double d = sum(pow3(weightedPSF(i, j, Range::all(), Range::all())));
                double b = c / d;
                weightedPSF(i, j, Range::all(), Range::all()) *= b;
                normalizationFactor(i, j) = d / (c * c);
            }
            
		}
    }
    
	biggerPSF.free();
        

}
//==============================================================================




//==============================================================================
/**
 * This function computes the weighted mask.
 * @param fluxMap Blitz 2-D array containing the flux of each pixel.
 */
void StepPhotometryOnTheGo::StepPhotometryOnTheGoMeasureStarsFluxes(Array<float, 2> fluxMap)
{
            
    subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY = p_DataSet->datasetGetsubFieldZeroY();
    
	//Initialization of the parameters to be measured and written into the photometry files (*.phot)
    magObs = 0.;
    fluxObs = 0.;
    magObsNorm = 0.;
    fluxObsNorm = 0.;
    magMC = 0.;
    fluxInput = 0.;
    magInput = 0.;
        

    //Auxiliar Blitz arrays initialization
	Array<float, 2> masks(fluxMap.shape());
	masks = 0.0;
	Array<float, 1> sortedPSF ;
	sortedPSF = 0.0;
	Array<int, 1> index1, index2;
    index1 = 0, index2 = 0;
    
	double realPSFCenter = (psfSubPixels - 1.) / (2. * subPixelsPerPixel);


	int x, y, counter, starCounter = 0, shiftPSFX, shiftPSFY;
	double backMedian = 0.0;
	double xd, yd, fractpartX, fractpartY, intpart;
	double psfSum = 0.0;
	double fluxSum = 0.0;
	
    //This threshold is to be used in the AP mask in such a way that this fraction 
    //of energy of the theoretical PSF lies within the aperture.

    double threshold = 0.9;    // AP photometry method takes 90% of the PSF energy.
    
    
	int cut = max(int(backgroundAnnulusOuterRadius), int(round(realPSFCenter)));
    
    //Check whether the BackgroundFlux must be automatically calculated as indicated in the XML 
    //input photometry parameters file. If so, photometryBackground < 0.
	if (photometryBackground < 0)
	{
		//determine median BACKGROUND in field
	    Array<float, 1> background;
		Array<float, 1> backArray(subPixelStarListOnCCD.rows());
		int backCounter = 0;
                    
        //For each star in the StarListOnCCD:
		for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
		{
            //Set x and y pixel coordinates of each star on the subpixel starList
			xd = (subPixelStarListOnCCD(i, 0) / subPixelsPerPixel) + subFieldZeroX;
			yd = (subPixelStarListOnCCD(i, 1) / subPixelsPerPixel) + subFieldZeroX;
    
			x = int(floor(xd));
			y = int(floor(yd));

            //If those coordinates are actually within the subfield frame:
			if (x > int(cut) && x < fluxMap.rows() - int(cut) && y > int(cut) && y < fluxMap.cols() - int(cut))
			{
				//compute background flux around aperture
				counter = 0;
                
                
                //Take a square field of backgroundAnnulusOuterRadius size
				background.resize(backgroundAnnulusOuterRadius*backgroundAnnulusOuterRadius);
                background = 0.0;
                
                //For every pixel in an squared frame of backgroundAnnulusOuterRadius size
				for (double a = -backgroundAnnulusOuterRadius; a <= backgroundAnnulusOuterRadius; a++)
                {
					for (double b = -backgroundAnnulusOuterRadius; b <= backgroundAnnulusOuterRadius; b++)
                    {
                        //If the pixel is between the outer ring and the inner ring:
						if (sqrt(a * a + b * b) >= backgroundAnnulusInnerRadius && sqrt(a * a + b * b) < backgroundAnnulusOuterRadius)
						{
                            //Add the flux of that pixel to the background 
							counter++;
							background(counter - 1) = double(fluxMap(x + int(a), y + int(b)));
						}
                    }
                }
				background.resizeAndPreserve(counter);
				backArray(backCounter) = MathTools::median(background);
                
				backCounter++;
			}
		}
        
		//we take as background value the median of all background measurements of all stars in the field. 
		backArray.resizeAndPreserve(backCounter);
        
		backMedian = MathTools::median(backArray);

        //Free arrays
		background.free();
		backArray.free();
        
	}
    
    //Else if, the background is provided in the XML photometry parameter file:
	else if (photometryBackground >= 0)
    {
        //Set the input background parameter.
		backMedian = photometryBackground;
    }
    
    
    
    //MEASURE THE FLUX
	double centralPixel = ceil(realPSFCenter);
	double distanceToCentralPixel = centralPixel - realPSFCenter;
	double fluxi;
    
	cut = int(round(realPSFCenter));
            
   
    //For each star encountered in the subfield, included in the subPixelStarListOnCCD array:
	for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
	{            
        //Calculation of the pixel coords of the current star from its subpixel coords (given at subPixelStarListOnCCD)
		xd = subPixelStarListOnCCD(i, 0) / subPixelsPerPixel;
		yd = subPixelStarListOnCCD(i, 1) / subPixelsPerPixel;
		
        //x, y star pixel coords
 		x = int(floor(xd));
		y = int(floor(yd));
 
        //Initialize starId parameter to identify each star
        int starId = subPixelStarListOnCCD(i, 3);
        
        //Check if the star is actually in the subfield, and not in the margin
		if (x > int(cut) && x < fluxMap.rows() - int(cut) && y > int(cut) && y < fluxMap.cols() - int(cut))
		{
			//compute the fractional part of the stellar coordinates to know how to shift the PSF
			fractpartX = modf(xd, &intpart) + distanceToCentralPixel;
			fractpartY = modf(yd, &intpart) + distanceToCentralPixel;
     
           
			//compute the number of sub-pixels to shift the PSF
			shiftPSFX = int(round(fractpartX * subPixelsPerPixel));
			shiftPSFY = int(round(fractpartY * subPixelsPerPixel));
             
            //Initialize parameters
            fluxSum = 0;  
            
            //Measure flux of each star in the subPixelStarListOnCCD array using Aperture method.
			if (photometryMethod == "AP")
			{
				//Determine which pixels of the rebinnedPSF are used for photometry for given threshold. Sort the rebinnedPSF,
				//then go from the highest value down until above threshold. Then add up all positions where the pixel value
				//is taken
				MathTools::shakersort(weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()), sortedPSF, index1, index2);
                
                //sortedPSF is sorted ascending
				int j = sortedPSF.size() - 1; 
				
				psfSum = 0;
                
                //compute stellar flux in aperture
				do
				{
					psfSum += sortedPSF(j);
					fluxi = fluxMap(x - int(centralPixel) + index1(j), y - int(centralPixel) + index2(j)) - backMedian;

                    //Check not to have negative values in the fluxMap map
					if (fluxi < 0.)
                    {
						fluxi = 0.;
                    }
                    
                    //Sum up all the fluxMap values 
					fluxSum += fluxi;
                    
//                  //The mask is not written to file in the photometry o-the-go as each mask is different  
//                    //For the fist exposure
////					if (frameIndex == 0)
//                    {   
//                        //Save the mask to be written to disk 
//                        masks(x - int(centralPixel) + index1(j), y - int(centralPixel) + index2(j)) = fluxi;
//                    }
                    
					j--;
				} while (j >= 0 && psfSum < threshold);
                
				//normalize fluxsum with the psfSum to get extrapolated flux (as if the complete PSF was used)
				normalizationFactor(shiftPSFX, shiftPSFY) = 1. / psfSum;
			}            

            
            
            //Measure flux of each star in the subPixelStarListOnCCD array using WEIGHTED MASK.
			else if (photometryMethod == "WM") 
			{
                //Measure flux of each star in the subPixelStarListOnCCD array using Weighted Mask algorithm as proposed by Reza Samadi.                                 
				fluxSum = sum(where(weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()) * 
                        (fluxMap(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), Range(y - int(centralPixel), y - int(centralPixel) +
                        weightedPSF.extent(fourthDim) - 1)) - backMedian) > 0., weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()) * 
                        (fluxMap(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), Range(y - int(centralPixel), y - int(centralPixel) + 
                        weightedPSF.extent(fourthDim) - 1)) - backMedian), 0.));
                            

//                  //The mask is not written to file in the photometry o-the-go as each mask is different                  
//                //For the fist exposure
//                if (frameIndex == 0)
//                {
//                    //Save the mask to be written to disk 
//					masks(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1)) +=
//                            (weightedPSF(shiftPSFX,shiftPSFY, Range::all(), Range::all()) * (fluxMap(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), 
//                            Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1)) - backMedian));
//                }			
            }
            
                        

	                
		//convert flux to magnitudes
			if (fluxSum > 0)
			{
				magObs = -2.5 * log10(fluxSum / fluxCorrectionFactor);
                
				magObsNorm = -2.5 * log10(fluxSum * normalizationFactor(shiftPSFX, shiftPSFY) / fluxCorrectionFactor);
                
				fluxObs = fluxSum * exposureTime;
				
                fluxObsNorm = fluxObs * normalizationFactor(shiftPSFX, shiftPSFY);
				
                magMC = Statistics::getPoisson(subPixelStarListOnCCD(i, 2), (unsigned int) time(0));
                
				magInput = subPixelStarListOnCCD(i, 2);
                
				fluxInput =  pow(10., -0.4 * magInput) * exposureTime * fluxCorrectionFactor;


				//write light curves to files
				stringstream fileStarOut;
                
				char tempChar[10];
				sprintf(tempChar, "%06d", starId);
				fileStarOut << photometryDirName + "/" << "/star" << tempChar << ".phot";
                
				ofstream outStar(fileStarOut.str().c_str(), std::ios_base::out | std::ios_base::app);
				outStar.precision(15);
                
                
                outStar << 0 * (exposureTime + readOutTime)<< " "; 

                outStar << magInput << " " << fluxInput;
                outStar << " " << magObs << " " << magObsNorm << " " << fluxObs << " " << fluxObsNorm << " " << backMedian 
                        * exposureTime << std::endl;
                    
  
				//write star coordinates to a file

				outStar.close();
                
				starCounter++;
                               
			}
		}
	}
    
    

	sortedPSF.free();
	index1.free();
	index2.free();
	masks.free();
    
}
//==============================================================================






//==============================================================================
/**
 * This function writes the photometry processing info to an ASCII file.
 * @param fileName File name to be written.
 */
void StepPhotometryOnTheGo::StepPhotometryOnTheGoWriteInfo(std::string fileName)
{
	ofstream out(fileName.c_str());
	if (!out.is_open())
	{
		std::cerr << "\nError (StepPhotometryOnTheGo::StepPhotometryOnTheGoWriteInfo()): Unable to open output information file " << fileName;
		exit(1);
	}
    
	std::string baseName = outputPath + "/" + prefix;
	baseName.erase(baseName.end() - 1);
    
	out << "-----------------------------------------------------------------------------" << std::endl;
	out << "Information about the photometry files in this directory" << std::endl;
	if (photometryMethod == "WM")   //If photometry method is Weighted Mask
		out << "PHOTOMETRY METHOD: WEIGHTED MASKS" << std::endl;
	else if (photometryMethod == "AP") //If photometry method is Aperture
		out << "PHOTOMETRY METHOD: APERTURE" << std::endl;
    
	if (flatfieldCorrection)
		out << "CORRECTED FOR FLATFIELD" << std::endl;
	else
		out << "NOT CORRECTED FOR FLATFIELD" << std::endl;
    
	if (frameTransferSmearingCorrection)
		out << "CORRECTED FOR TRAILING" << std::endl;
	else
		out << "NOT CORRECTED FOR TRAILING" << std::endl;
    
	out << "Number of Telescopes (N): " << photometryNumTelescopes << std::endl;
	out << "-----------------------------------------------------------------------------" << std::endl;
	out << baseName + ".info" << "  ............. Input parameters for simulation" << std::endl;
	out << baseName + "*.fits" << "  ............. Raw sub-images" << std::endl;
	out << baseName + "_flat.fits" << "  ............. Flatfield (illumination flux 30000 e)" << std::endl;
	out << baseName + "_reduced0.fits" << "  ............. First image corrected by trailing, flatfield, gain and exposure (flux/s)"
    << std::endl;
	out << baseName + "_masks0.fits" << " ................ Masks of the first image" << std::endl;
	out << "Normalized flux is computed from the theoretical flux of the model PSF at the sub-pixel position of the star" << std::endl;

	out << "star*.phot" << " ............. Lightcurve of each star" << std::endl;
	out << "                                      Column description:" << std::endl;
	out << "                                      1: time [s]" << std::endl;
	out << "                                      2: Input magnitude" << std::endl;
	out << "                                      3: Input flux" << std::endl;
	out << "                                      4: Measured magnitude" << std::endl;
	out << "                                      5: Measured normalized magnitude" << std::endl;
	out << "                                      6: Measured flux" << std::endl;
	out << "                                      7: Measured normalized flux" << std::endl;
	out << "                                      8: Background flux" << std::endl;
    
	out << "star_coordinates.dat" << " ............. Coordinates and magnitude of measured stars" << std::endl;
	out << "                                      Column description:" << std::endl;
	out << "                                      1: File name with lightcurve" << std::endl;
	out << "                                      2: X-coordinate of star on sub-image" << std::endl;
	out << "                                      3: Y-coordinate of star on sub-image" << std::endl;
	out << "                                      4: Input magnitude" << std::endl;
	out << "                                      5: Input flux" << std::endl;
    
	
    
	out << baseName + "_psf.dat" << " ............. Input PSF of image" << std::endl;
	out << "                                      Column description:" << std::endl;
	out << "                                      1: X-coordinate" << std::endl;
	out << "                                      2: Y-coordinate" << std::endl;
	out << "                                      3: Normalized flux" << std::endl;
    
	
	out.close();
}
//==============================================================================





