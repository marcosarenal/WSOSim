///////////////////////////////////////////////////////////
//  ProcessingPhotometry.cpp
//  Implementation of the Class ProcessingPhotometry
//  Created on:      23-Oct-2012 2:00:00 PM
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



#include "ProcessingPhotometry.h"


//==============================================================================
/**
 * Constructor method
 */
ProcessingPhotometry::ProcessingPhotometry(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ProcessingPhotometry::~ProcessingPhotometry()
{
    normalizationFactor.free();
}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function performs the Photometry processing reading the FITS files. 
 * The flux of each star in the 
 * sub-field is measured assuming a Gaussian weighted mask and then derived the 
 * noise-to-signal ratio (NSR).
 */
void ProcessingPhotometry::processingPhotometryPipeline(DataSet &m_DataSet, DataSetPhotometry &datasetPhotometry)
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
        
    //Retrieving parameters from PhotometryDataSet.
    photometryDirName = datasetPhotometry.datasetGetphotometryDirName();
    photometryPlotsDir = datasetPhotometry.datasetGetphotometryPlotsDir();
    photometryNumTelescopes = datasetPhotometry.datasetGetphotometryNumTelescopes();
   
    photometryMethod = datasetPhotometry.datasetGetphotometryMethod();
    photometryBackground = datasetPhotometry.datasetGetphotometryBackground();
    
    frameTransferSmearingCorrection = datasetPhotometry.datasetGetframeTransferSmearingCorrection();
    flatfieldCorrection = datasetPhotometry.datasetGetphotometryFlatfieldCorrectionr();
    
    backgroundAnnulusInnerRadius = datasetPhotometry.datasetGetphotometryBackgroundAnnulusInnerRadius();
    backgroundAnnulusOuterRadius = datasetPhotometry.datasetGetphotometryBackgroundAnnulusOuterRadius();
       
    
    //Retrieving the subPixelStarListOnCCD contains for each star, its X and Y position in CCD, magnitude, RA, 
    //declination and identification number.
    subPixelStarListOnCCD.resize(m_DataSet.datasetGetsubPixelStarListOnCCD().extent(0),m_DataSet.datasetGetsubPixelStarListOnCCD().extent(1));
    subPixelStarListOnCCD = 0.0;
    subPixelStarListOnCCD = m_DataSet.datasetGetsubPixelStarListOnCCD();

    
    
    //Recompute the star coordinates on the CCD to take the cutting of the edges into account.
    ProcessingPhotometry::processingPhotometrySubPixelStarListOnCCDCut();
    
   
        
    //Compute the Weighted PSF mask.
    ProcessingPhotometry::processingPhotometryPreComputeWeightedMask();
    
                
    LogManager::log <<"    Number of frames: " << exposuresNamesArray.size();
    GlobalVariables::logManager.LogManagerShowLog();    

	//Define a progress counter variable
    double percent;
    percent = 0.0;
    
    //Iteration for each FITS image
    for (unsigned int iterExposure = 0; iterExposure < exposuresNamesArray.size(); iterExposure++)
    {
        //Just a progress counter
        percent = 100. * iterExposure / double(exposuresNamesArray.size());
        double dummy1 = fmod(percent, 10);
        double dummy2 = fmod(100. * (iterExposure - 1) / double(exposuresNamesArray.size()), 10);

        if (dummy1 <= dummy2)
        {
            LogManager::log << "    "<<percent << "%  processed";
            GlobalVariables::logManager.LogManagerShowLog();
        }

        //Take the corresponding FITS file
        string fitsFile = outputPath + "/" + prefix + "/" +  prefix + exposuresNamesArray(iterExposure) + ".fits";

        //Check that it is correctly opened
        if (FileUtilities::fileExists(fitsFile))
            {                       
            //Read the science image and set it in the rawImage array
            FileUtilities::FileUtilitiesReadFits(outputPath + "/" + prefix + "/" +  prefix + exposuresNamesArray(iterExposure) + ".fits", rawImage); 

            //Cut the pre and overscans off the science image
            flux.resize(rawImage.rows(), rawImage.cols() - numPrescanRows - numSmearingOverscanRows);
            flux = 0.0;

            //Set the flux values as taken from the science FITS files.
            flux = rawImage(Range::all(), Range(numPrescanRows, rawImage.cols() - numSmearingOverscanRows - 1));


            //Obtain the star smearing (trailing) register from the over-scan
            Array<float, 2> trailingMap;
            trailingMap.resize(rawImage.rows(), numSmearingOverscanRows);
            trailingMap = 0.0;
            trailingMap = rawImage(Range::all(), Range(rawImage.cols() - numSmearingOverscanRows, rawImage.cols() - 1));


            //TODO: correct BIAS from Prescan values instead of directly from the BIAS included in the overscan (pma))
            //Obtain the bias register from the pre-scan
            //            Array<float, 2> biasMap;
            //			biasMap.resize(rawImage.rows(), numPrescanRows);
            //            biasMap = 0.0;
            //			biasMap = rawImage(Range::all(), Range(0,  numPrescanRows));


            //Remove trailing from flux map
            if (frameTransferSmearingCorrection && numSmearingOverscanRows > 0)
            {
                for (int j = 0; j < trailingMap.rows(); j++)
                {    
                    flux(j, Range::all()) -= MathTools::median(trailingMap(j, Range::all()));            
                }
            }

            //After subtraction of trails, flux might be negative! Set all negative values to 0
            flux = where(flux < 0., 0., flux);


            //correct with flatfield
            //				if (correctFlatfield)
            //                {
            //m_DataSet.datasetGetFlatFieldMap(flatfieldMap);
            flatfieldMap.resize(flux.rows(), flux.cols());
            flatfieldMap = 0.0;
            
            flatfieldMap = m_DataSet.datasetGetFlatFieldMap();

            flux /= flatfieldMap;
            //                }

            //Correct flux with gain, exposure time 
            flux *= (gain / (exposureTime ));

            //Write the corrected flux map and the trailing map to FITS file only for the first exposure
            if (iterExposure == 0)
            {              
                FileUtilities::FileUtilitiesWriteFits(photometryDirName + "/" + prefix + "_reduced0", flux);
                FileUtilities::FileUtilitiesWriteFits(photometryDirName + "/" + prefix + "_trailing0", trailingMap, 0);                    
            }


            //Compute the photometry mask
            ProcessingPhotometry::processingPhotometryMeasureStarsFluxes(flux, iterExposure);	

            //Free the created Map
            trailingMap.free();
            //            biasMap.free();


        }
    }

    ProcessingPhotometry::processingPhotometryComputeStatistics(photometryPlotsDir + "/" + prefix  + ".stat");

    //write psf as ascii file
    ProcessingPhotometry::processingPhotometryWritePsfToASCIIFile(photometryPlotsDir + "/" + prefix  + "_psf.dat");
    ProcessingPhotometry::processingPhotometryWriteInfo(photometryPlotsDir + "/README");
    
}
//==============================================================================







//==============================================================================
/**
 * This function takes the subPixelStarListOnCCD array and recompute the star coordinates on the CCD to take the cutting of the edges into account.
 */
void ProcessingPhotometry::processingPhotometrySubPixelStarListOnCCDCut()
{
    Array<float, 2>  tempSubPixelStarListOnCCD;
    
    tempSubPixelStarListOnCCD.resize(subPixelStarListOnCCD.rows(), subPixelStarListOnCCD.cols());
    tempSubPixelStarListOnCCD = 0.0;

	int starCounter = 0;
	int x, y;
        
	for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
	{
		x = int(round(subPixelStarListOnCCD(i, 0)));
		y = int(round(subPixelStarListOnCCD(i, 1)));
		if (x >= edgePixels * subPixelsPerPixel  && x < (subFieldSizeX - edgePixels) * subPixelsPerPixel  && 
                    y >= edgePixels * subPixelsPerPixel && y < (subFieldSizeY - edgePixels) * subPixelsPerPixel)
		{
			tempSubPixelStarListOnCCD(starCounter, 0) = subPixelStarListOnCCD(i, 0) - double(edgePixels * subPixelsPerPixel);
			tempSubPixelStarListOnCCD(starCounter, 1) = subPixelStarListOnCCD(i, 1) - double(edgePixels * subPixelsPerPixel);
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
    
    //Free array
    tempSubPixelStarListOnCCD.free();
}
//==============================================================================




//==============================================================================
/**
 * This function computes the weighted mask using the algorithm as proposed by Reza Samadi.
 */
void ProcessingPhotometry::processingPhotometryPreComputeWeightedMask()
{
    
    //Retrieving parameters from DataSet
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    psfNumPixels = p_DataSet->datasetGetpsfNumPixels();


    //Create a bigger PSF map
	Array<float, 2> biggerPSF(psfMap.rows() + 2 * subPixelsPerPixel -1 , psfMap.rows() + 2 * subPixelsPerPixel -1);
    
    //Initialize the weighted mask
	weightedPSF.resize(subPixelsPerPixel + 1, subPixelsPerPixel + 1, int(biggerPSF.rows() / subPixelsPerPixel), int(biggerPSF.cols() / subPixelsPerPixel));
    weightedPSF = 0.0;
           
    //Resize the normalization factor matrix
	normalizationFactor.resize(subPixelsPerPixel + 1, subPixelsPerPixel + 1);
    normalizationFactor= 0.0;
 
    //For each subpixel in the weighted PSF mask
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
         
    //Free the generated array
	biggerPSF.free();
        

}
//==============================================================================








//==============================================================================
/**
 * This function computes the weighted mask.
 * @param fluxMap Blitz 2-D array containing the flux of each pixel.
 * @param frameIndex Number of frame to be processed.
 */
void ProcessingPhotometry::processingPhotometryMeasureStarsFluxes(Array<float, 2> flux, int frameIndex)
{
	//Initialization of the parameters to be measured and written into the photometry files (*.phot)
    magObs = 0.;
    fluxObs = 0.;
    magObsNorm = 0.;
    fluxObsNorm = 0.;
    magMC = 0.;
    fluxInput = 0.;
    magInput = 0.;
                    
    //Auxiliar Blitz arrays initialization
	Array<float, 2> masks(flux.shape());
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
                

    Array<double, 1> sortedSmallMask ;
    sortedSmallMask = 0.0;
                
                
    //Define the file to store the star coordinates data
    string fileCoordOutString(photometryDirName + "/star_coordinates.dat");
    ofstream outCoorFile(fileCoordOutString.c_str(), std::ios_base::out | std::ios_base::app);
    outCoorFile.precision(15);

                    

	int cut = max(int(backgroundAnnulusOuterRadius), int(round(realPSFCenter)));

    //Check whether the BackgroundFlux must be automatically calculated as indicated in the XML 
    //input photometry parameters file. If so, photometryBackground < 0.
	if (photometryBackground < 0)
	{
		//BACKGROUND CALCULATION
        //Take a square field of backgroundAnnulusOuterRadius size
	    Array<float, 1> background(backgroundAnnulusOuterRadius*backgroundAnnulusOuterRadius);
        background = 0.0;
		Array<float, 1> backArray(subPixelStarListOnCCD.rows());
		int backCounter = 0;
        int backsize = 1000;
                     
        //For each star in the StarListOnCCD:
		for (int i = 0; i < subPixelStarListOnCCD.rows(); i++)
		{
            //Set x and y pixel coordinates of each star on the subpixel starList
			xd = subPixelStarListOnCCD(i, 0) / subPixelsPerPixel;
			yd = subPixelStarListOnCCD(i, 1) / subPixelsPerPixel;
    
			x = int(floor(xd));
			y = int(floor(yd));
                          
            //If those coordinates are actually within the subfield frame:
			if (x > int(cut) && x < flux.rows() - int(cut) && y > int(cut) && y < flux.cols() - int(cut))
			{  
				//compute background flux around aperture
				counter = 0;
                
                //resize background map
                background.resize(backsize + 100);
                           
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
							background(counter - 1) = double(flux(x + int(a), y + int(b)));
                                    
						}
                    }
            
                }
				background.resizeAndPreserve(counter);
				backsize = counter;
				backArray(backCounter) = MathTools::median(background);
                                                    
				backCounter++;
			}
		}
     
		background.free();
        
        //we take as background value the median of all background measurements of all stars in the field
		backArray.resizeAndPreserve(backCounter);

        //Check whether there is background flux detected
        if (backArray.rows() != 0)
        {
            backMedian = MathTools::median(backArray);
            if (backMedian == 0)
            {
                cerr << "\nError (ProcessingPhotometry::processingPhotometryMeasureStarsFluxes()): No background detected on images."<<endl;
                exit (1);
            }
        }
        else
        {      
            cerr << "\nError (ProcessingPhotometry::processingPhotometryMeasureStarsFluxes()): No background detected on images."<<endl;
            exit (1);
        }
        
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
    //END background calculation                                                                          


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
        
        //Check if the star is actually in the subfield, and not in the margin
		if (x > int(cut) && x < flux.rows() - int(cut) && y > int(cut) && y < flux.cols() - int(cut))
		{
			//compute the fractional part of the stellar coordinates to know how to shift the PSF
			fractpartX = modf(xd, &intpart) + distanceToCentralPixel;
			fractpartY = modf(yd, &intpart) + distanceToCentralPixel;
     
           
			//compute the number of sub-pixels to shift the PSF
			shiftPSFX = int(round(fractpartX * subPixelsPerPixel));
			shiftPSFY = int(round(fractpartY * subPixelsPerPixel));
			
            //Initialize parameter
            fluxSum = 0.0;

            //Measure flux of each star in the subPixelStarListOnCCD array using Aperture method.
			if (photometryMethod == "AP")
			{                
				//Determine which pixels of the rebinned PSF are used for photometry for given threshold. Sort the rebinned PSF,
				//then go from the highest value down until above threshold. Then add up all positions where the pixel value
				//is taken.
				MathTools::shakersort(weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()), sortedPSF, index1, index2);

                //sortedPSF is sorted ascending
				int j = sortedPSF.size() - 1; 
                
				psfSum = 0.0;

				//compute stellar flux in aperture
				do
				{
					psfSum += sortedPSF(j);
					fluxi = flux(x - int(centralPixel) + index1(j), y - int(centralPixel) + index2(j)) - backMedian;
					
                    //Check not to have negative values in the flux map
                    if (fluxi < 0.)
                    {
						fluxi = 0.;
                    }
                    
                    //Sum up all the flux values 
					fluxSum += fluxi;

                    //For the fist exposure
					if (frameIndex == 0)
                    {
                        //Save the mask to be written to disk 
						masks(x - int(centralPixel) + index1(j), y - int(centralPixel) + index2(j)) = fluxi;
                    }
					j--;
				} while (j >= 0 && psfSum < threshold);
                
				//normalize fluxsum with the psfSum to get extrapolated flux (as if the complete PSF was used)
				normalizationFactor(shiftPSFX, shiftPSFY) = 1. / psfSum;   


			}            
            
            //Measure flux of each star in the subPixelStarListOnCCD array using WEIGHTED MASK.
			else if (photometryMethod == "WM") 
			{      
                //Measure flux of each star in the subPixelStarListOnCCD array using Weighted Mask algorithm as proposed by Reza Samadi.
                //Note:where(array-expr1, array-expr2, array-expr3)
                //     Wherever array-expr1 is true, array-expr2 is returned. Where array-expr1 is false, array-expr3 is returned. 
  				fluxSum = sum(where(
                                weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()) * 
                                (   flux(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), 
                                            Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1))
                                    - backMedian) > 0., 
                                weightedPSF(shiftPSFX, shiftPSFY, Range::all(), Range::all()) * 
                                (   flux(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), 
                                            Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1))
                                    - backMedian), 
                                0.0)
                             ); 
                                         
                //For the fist exposure
                if (frameIndex == 0)
                {
                    //Save the mask to be written to disk 
					masks(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), 
                          Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1)) +=
                            (weightedPSF(shiftPSFX,shiftPSFY, Range::all(), Range::all()) * 
                            (flux(Range(x - int(centralPixel), x - int(centralPixel) + weightedPSF.extent(thirdDim) - 1), 
                                     Range(y - int(centralPixel), y - int(centralPixel) + weightedPSF.extent(fourthDim) - 1)) - backMedian));
                }
            }
            //If none of the available masks is selected
            else
            {                
                cerr << "\nError (ProcessingPhotometry::processingPhotometryMeasureStarsFluxes()): "<< photometryMethod<<"  is not a valid photometry method."<<endl;
                exit (1);
                
            }        
	                
		    //convert flux to magnitudes
			if (fluxSum > 0)
			{            
                                            
				magObs = -2.5 * log10(fluxSum / (fluxm0 * areaTelescope * transEff * quantEff));

                
				magObsNorm = -2.5 * log10(fluxSum * normalizationFactor(shiftPSFX, shiftPSFY) / (fluxm0 * areaTelescope * transEff * quantEff));
                
				fluxObs = fluxSum * exposureTime;
				
                fluxObsNorm = fluxObs * normalizationFactor(shiftPSFX, shiftPSFY);
				
                magMC = -2.5 * log10(Statistics::getPoisson(subPixelStarListOnCCD(i, 3), (unsigned int) time(0))
                                     / (exposureTime * fluxm0 * areaTelescope * transEff * quantEff));
                
				magInput = subPixelStarListOnCCD(i, 2);
                
				fluxInput = fluxm0 * exposureTime * pow(10., -0.4 * magInput) * areaTelescope * transEff * quantEff;

                               
                
				//write light curves to files
				stringstream fileStarOut;
                
				char tempChar[10];
				sprintf(tempChar, "%06d", starCounter);
				fileStarOut << photometryDirName + "/" << "/star" << tempChar << ".phot";
                
				ofstream outStar(fileStarOut.str().c_str(), std::ios_base::out | std::ios_base::app);
				outStar.precision(15);
                
                stringstream ss;
                ss << photometryDirName << "/star_input" << i << ".dat";
                
                string inputMagFile = ss.str();
                                
                outStar << frameIndex * (exposureTime + readOutTime)<< " "; 
                if (FileUtilities::fileExists(inputMagFile))
                {
                    outStar << FileUtilities::readValue(inputMagFile, 3, frameIndex) << " " << FileUtilities::readValue(inputMagFile, 4, frameIndex);
                }
                else
                {
                    outStar << magInput << " " << fluxInput;
                    outStar << " " << magObs << " " << magObsNorm << " " << fluxObs << " " << fluxObsNorm << " " << backMedian 
                            * exposureTime << endl;
                    
                }
                
				//write star coordinates to a file
				if (frameIndex == 0)
				{                
					outCoorFile << "star" << starCounter << ".phot" << " " << xd << " " << yd;
					outCoorFile << " " << subPixelStarListOnCCD(i, 2) << " " << subPixelStarListOnCCD(i, 3)
                    << endl;
				}                				
                outStar.close();

				starCounter++;
                                               
			}				
		}					
    
	}    
    //Write to FITS the Photometry mask of the first frame
	if (frameIndex == 0)
	{
        FileUtilities::FileUtilitiesWriteFits(photometryDirName + "/" + prefix + "_masks0", masks, 0);

	}        
                   
    //Close the photometry file 
    outCoorFile.close();                
   
                
}
//==============================================================================





//==============================================================================
void ProcessingPhotometry::processingPhotometryComputeStatistics(string fileName)
{
    LogManager::log << "    Computing statistics..." ;
    GlobalVariables::logManager.LogManagerShowLog();
    
	ofstream out(fileName.c_str());
	if (!out.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryComputeStatistics()): Unable to open output information file";
		exit(1);
	}
	out.precision(15);
    
	out << "### File Name | Input mag | Input F | Mean input mag (Monte-Carlo) | Std. dev of input mag (Monte-Carlo) |";
	out << " Mean measured mag | Mean measured normalized mag | Std. dev of measured mag | RMS of measured mag | Mean measured F |";
	out << " Mean measured normalized flux | Std. dev of measured flux | Theoretical photon noise (flux) | Measured noise ppm/hour for "
    << photometryNumTelescopes << " telescopes |";
	out << " Theoretical photon noise ppm/hour for " << photometryNumTelescopes << " telescopes | Measured N (8T) | Theoretical PN (8T) |";
	out << " Measured N (16T) | Theoretical PN (16T) | Measured N ("<<photometryNumTelescopes<<") | Theoretical PN ("<<photometryNumTelescopes<<")" << endl;
	double stdDevMag, medianMag, stdDevInputMag, medianInputMag, rmsMag, theoreticalPhotonNoise, stdDevFlux, medianFlux, ppmphObs, ppmphPN,
    medianFluxNorm, medianMagNorm;
    
	Array<float, 1> magObsArray, fluxObsArray, magObsNormArray, fluxObsNormArray;
	Array<float, 1> magMCArray, fluxInputArray;
	double medianFluxInput, medianMagInput;
    
	//compute number of stars(==number of files containing star*.phot)
	vector<string> starFiles;
	FileUtilities::getDir(photometryDirName, starFiles, ".phot");
	int numStars = starFiles.size();
	MathTools::shakersort(starFiles);
    
	double normFactor;
	string starPath;
    
    LogManager::log << "    Number of stars on frame: " << numStars ;
    GlobalVariables::logManager.LogManagerShowLog();
    
	for (int i = 0; i < numStars; i++)
	{
		starPath = photometryDirName + "/" + starFiles[i];
		FileUtilities::readCol(starPath, 1, magMCArray, 0);
		FileUtilities::readCol(starPath, 2, fluxInputArray, 0);
		FileUtilities::readCol(starPath, 3, magObsArray, 0);
		FileUtilities::readCol(starPath, 4, magObsNormArray, 0);
		FileUtilities::readCol(starPath, 5, fluxObsArray, 0);
		FileUtilities::readCol(starPath, 6, fluxObsNormArray, 0);
		medianFluxInput = MathTools::median(fluxInputArray);
		medianMagInput = MathTools::median(magMCArray);
        
		stdDevMag = Statistics::getStdDev(magObsArray); //compute the standard deviation of the measured magnitude

		medianMag = MathTools::median(magObsArray);        
        
		medianMagNorm = MathTools::median(magObsNormArray); //take the median of the normalized magnitude to be able to compare with input magnitude
		stdDevFlux = Statistics::getStdDev(fluxObsArray);
		medianFlux = MathTools::median(fluxObsArray);
		medianFluxNorm = MathTools::median(fluxObsNormArray);
		normFactor = medianFluxNorm / medianFlux;
		stdDevInputMag = Statistics::getStdDev(magMCArray);
		medianInputMag = MathTools::median(magMCArray);
		rmsMag = sqrt(pow2(medianMag) + pow2(stdDevMag));
		theoreticalPhotonNoise = sqrt(medianFlux);
		ppmphObs = (stdDevFlux * normFactor / medianFluxInput) * 1.0e6 / sqrt(3600. / exposureTime); //measured noise ppm per hour (stdDevFlux/medianFlux)
		ppmphPN = (1. / theoreticalPhotonNoise) * 1.0e6 / sqrt(3600. / exposureTime); //theoretical photon noise ppm per hour
		out << starFiles[i] << " ";
		out << medianMagInput << " " << medianFluxInput << " " << medianInputMag;
		out << " " << stdDevInputMag << " " << medianMag << " " << medianMagNorm << " " << stdDevMag << " " << rmsMag << " ";
		out << medianFlux << " " << medianFluxNorm << " " << stdDevFlux << " " << theoreticalPhotonNoise << " ";
		out << ppmphObs << " " << ppmphPN << " " << ppmphObs / sqrt(photometryNumTelescopes) << " " << ppmphPN / sqrt(photometryNumTelescopes) << " " << ppmphObs / sqrt(8)
        << " " << ppmphPN / sqrt(8) << " ";
		out << ppmphObs / sqrt(16) << " " << ppmphPN / sqrt(16) << " " << ppmphObs / sqrt(24) << " " << ppmphPN / sqrt(24) << " "
        << ppmphObs / sqrt(photometryNumTelescopes) << " " << ppmphPN / sqrt(photometryNumTelescopes) << endl;
	}
    
	ProcessingPhotometry::processingPhotometryMakeAnalysis(fileName);
    
	magObsArray.free();
	fluxObsArray.free();
	magObsNormArray.free();
	fluxObsNormArray.free();
	magMCArray.free();
	fluxInputArray.free();
}
//==============================================================================




//==============================================================================
void ProcessingPhotometry::processingPhotometryMakeAnalysis(string fileName)
{

	Array<float, 1> magInput, magMeasured, fluxInput;
	FileUtilities::readCol(fileName, 1, magInput, 0);
	FileUtilities::readCol(fileName, 6, magMeasured, 0);
	FileUtilities::readCol(fileName, 2, fluxInput, 0);

	Array<float, 1> magDiff(magInput.shape());
	magDiff = magMeasured - magInput;

	string fileOut = photometryPlotsDir + "/crowding_" + prefix + +".dat";
    
	ofstream out(fileOut.c_str());
	if (!out.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryMakeAnalysis()): Unable to open " << fileOut << endl;
		exit(1);
	}
	out.precision(15);
    
	double counter1, counter2, ratio;
    
	for (double f = 3; f < 18; f += 1)
	{
		counter1 = 0;
		counter2 = 0;
		ratio = 0.;
		for (unsigned int i = 0; i < magInput.size(); i++)
		{
			if (magInput(i) >= f && magInput(i) < f + 1.)
				counter1++;
			if (magInput(i) >= f && magInput(i) < f + 1. && magDiff(i) < -0.113) // 1/3 noise criterion
				counter2++;
		}
		if (counter1 > 0)
			ratio = 100. - 100. * counter2 / counter1;
		out << f + 0.5 << " " << 1. << " " << magInput.size() << " " << counter1 << " " << counter2 << " " << ratio << endl;
	}
    
    
	//Compute mean/median noise at each magnitude
    
	Array<float, 1> noise;
	FileUtilities::readCol(fileName, 13, noise, 0);
    
	string noiseFileName = photometryPlotsDir + "/noise_" + prefix + ".dat";
    
	ofstream out2(noiseFileName.c_str());
	if (!out2.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryMakeAnalysis()): Unable to open "  << noiseFileName << endl;
		exit(1);
	}
	out2.precision(15);
    
	out2 << "### Magnitude | Size of bin | Total no. of stars | No. of stars in bin | Mean noise (T=1) | Median (T=1) | Mean noise (T="
    << photometryNumTelescopes << ") | Median (T=" << photometryNumTelescopes << ") ";
	out2
    << "| Mean noise (T=8) | Median (T=8) | Mean noise (T=16) | Median (T=16) | Mean noise (T=24) | Median (T=24)| Mean noise (T="<<photometryNumTelescopes<<") | Median (T="<<photometryNumTelescopes<<")"
    << endl;
    
	Array<float, 1> noisePerMag;
	double meany, mediany, dmag;
	dmag = 0.5;
    
	for (double f = 3 + dmag / 2.; f < 18; f += dmag)
	{
		noisePerMag.resize(100000);
		counter1 = 0;
		ratio = 0.;
		for (unsigned int i = 0; i < magInput.size(); i++)
		{
			if (magInput(i) >= f && magInput(i) < f + dmag)
			{
				noisePerMag(counter1) = noise(i);
				counter1++;
			}
		}
		noisePerMag.resizeAndPreserve(counter1);
		if (counter1 > 0)
		{
			meany = mean(noisePerMag);
			mediany = MathTools::median(noisePerMag);
			out2 << f + dmag / 2. << " " << dmag << " " << magInput.size() << " " << counter1 << " " << meany << " " << mediany;
			out2 << " " << meany / sqrt(photometryNumTelescopes) << " " << mediany / sqrt(photometryNumTelescopes) << " " << meany / sqrt(8.) << " " << mediany
            / sqrt(8.) << " " << meany / sqrt(16.) << " " << mediany / sqrt(16.);
			out2 << " " << meany / sqrt(24.) << " " << mediany / sqrt(24.) << " " << meany / sqrt(photometryNumTelescopes) << " " << mediany / sqrt(photometryNumTelescopes)
            << endl;
		}
	}
    
	//Compute number of stars below certain threshold levels
	Array<int, 1> counter_ppmX(20);
    
	string fileOut3;
	stringstream sstr;
	sstr << photometryPlotsDir << "/starcounts_" << photometryNumTelescopes << "T_" << prefix << ".dat";
	fileOut3 = sstr.str();
    
	ofstream out3(fileOut3.c_str());
	if (!out3.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryMakeAnalysis()): Unable to open "  << fileOut3 << endl;
		exit(1);
	}
	out3.precision(15);
    
	out3 << "### Noise level  |  N(stars) <= noise level ppm/h | % stars <= noise level ppm/h" << endl;

	int counter_ppm27 = 0;
	int counter_ppm80 = 0;
    
	counter_ppmX = 0;
    
	noise /= sqrt(photometryNumTelescopes);
    
	for (unsigned int i = 0; i < noise.size(); i++)
	{
		if (noise(i) <= 27)
			counter_ppm27++;
		if (noise(i) <= 80)
			counter_ppm80++;
		for (unsigned int j = 1; j < counter_ppmX.size(); j++)
			if (noise(i) <= j * 10)
				counter_ppmX(j)++;
	}
	out3 << 27 << " " << counter_ppm27 << " " << 100. * counter_ppm27 / noise.size() << endl;
	out3 << 80 << " " << counter_ppm80 << " " << 100. * counter_ppm80 / noise.size() << endl;
	for (unsigned int j = 1; j < counter_ppmX.size(); j++)
		out3 << j * 10 << " " << counter_ppmX(j) << " " << 100. * counter_ppmX(j) / noise.size() << endl;
    
	//++++++++++++++++++++++++++++++++++++++++
	string gnu = photometryPlotsDir + "/gnuscript.gnu";
	ofstream outGnu(gnu.c_str());
	if (!outGnu.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryMakeAnalysis()): output gnuscript file " << gnu;
		exit(1);
	}
    
    LogManager::log << "    Creating script for gnuplot (" << gnu << ")" ;
    GlobalVariables::logManager.LogManagerShowLog();
    
	outGnu.precision(15);
	outGnu << "#!/bin/sh" << endl;
	outGnu << "gnuplot << EOF" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "crowding.jpg'" << endl;
	outGnu << "set key bottom" << endl;
	outGnu << "f(x)=x" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "set ylabel 'Measured magnitude'" << endl;
	outGnu << "set style fill solid border -1" << endl;
	outGnu << "plot '" << fileName << "' u 2:7 title '', f(x)" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "crowding.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "crowding_stat.jpg'" << endl;
	outGnu << "set xlabel 'Magnitude'" << endl;
	outGnu << "set ylabel 'Occurence [%]'" << endl;
	outGnu << "set xrange [5:17]" << endl;
	outGnu << "plot '" << fileOut << "' u 1:6 w boxes title ''" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "crowding_stat.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	/* outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	 outGnu << "set output '" << outputPath + "/" + FileUtilities::getCompleteBaseName(fileName) << "variables_stat.jpg'" << endl;
	 outGnu << "set xlabel 'Magnitude'" << endl;
	 outGnu << "set ylabel 'Occurence [%]'" << endl;
	 outGnu << "set key top" << endl;
	 outGnu << "set boxwidth 1"<<endl;
	 outGnu << "plot '" << fileOut2.fileName() << "' u 1:10 w boxes title '10ppm puls. detectable (single star puls)','";
	 outGnu << fileOut2.fileName() << "' u 1:11:2 w boxes title '10ppm puls. detectable (50% stars puls)','";
	 outGnu << fileOut2.fileName() << "' u 1:12:3 w boxes title '10ppm puls. detectable (100% stars puls)'" << endl;
	 outGnu << "set terminal post landscape color solid enhanced 13" << endl;
	 outGnu << "set output '" << outputPath + "/" + FileUtilities::getCompleteBaseName(fileName) << "variables_stat.eps'" << endl;
	 outGnu << "rep" << endl;
	 outGnu << "set out" << endl;
	 */
    
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmph.jpg'" << endl;
	outGnu << "set key bottom" << endl;
	outGnu << "set log y" << endl;
	outGnu << "set yrange [1:]" << endl;
	outGnu << "set xrange [*:*]" << endl;
	outGnu << "f(x)=27." << endl;
	outGnu << "g(x)=80.0" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Measured magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 6:14 title 'observed noise', '" << fileName
    << "' u 6:15 w l title 'theoretical photon noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmph.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag" << photometryNumTelescopes << "tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:(\\$14/sqrt(" << photometryNumTelescopes << ")) title 'observed noise', '" << fileName
    << "' u 6:(\\$15/sqrt(" << photometryNumTelescopes << ")) w l title 'theoretical photon noise','" << noiseFileName
    << "' u 2:9 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag" << photometryNumTelescopes << "tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag1tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:14 title 'observed noise', '" << fileName << "' u 6:15 w l title 'theoretical photon noise','"
    << noiseFileName << "' u 2:7 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag1tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag8tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:18 title 'observed noise', '" << fileName << "' u 6:19 w l title 'theoretical photon noise','"
    << noiseFileName << "' u 1:10 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag8tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag16tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:20 title 'observed noise', '" << fileName << "' u 6:21 w l title 'theoretical photon noise','"
    << noiseFileName << "' u 1:12 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag16tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag24tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:22 title 'observed noise', '" << fileName << "' u 6:23 w l title 'theoretical photon noise','"
    << noiseFileName << "' u 1:14 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag24tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag"<<photometryNumTelescopes<<"tel.jpg'" << endl;
	outGnu << "set ylabel 'ppm/hour'" << endl;
	outGnu << "set xlabel 'Input magnitude'" << endl;
	outGnu << "plot '" << fileName << "' u 2:24 title 'observed noise', '" << fileName << "' u 6:25 w l title 'theoretical photon noise','"
    << noiseFileName << "' u 1:16 w l title 'median obs. noise', f(x) title '27 ppm/hour', g(x) title '80.0 ppm/hour'" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "ppmphinputmag"<<photometryNumTelescopes<<"tel.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
	outGnu << "reset" << endl;
	outGnu << "set terminal jpeg large size 1000,800 enhanced" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "psf.jpg'" << endl;
	outGnu << "set ylabel 'y'" << endl;
	outGnu << "set xlabel 'x'" << endl;
	outGnu << "set pm3d map" << endl;
	outGnu << "set xrange[-4:4]" << endl;
	outGnu << "set yrange[-4:4]" << endl;
	outGnu << "splot '" << photometryPlotsDir + "/" + prefix + "_psf.dat" << "' title ''" << endl;
	outGnu << "set terminal post landscape color solid enhanced 15" << endl;
	outGnu << "set output '" << photometryPlotsDir + "/" + prefix << "psf.eps'" << endl;
	outGnu << "rep" << endl;
	outGnu << "set out" << endl;
    
	outGnu << "set nomultiplot\nset terminal x11\nset output" << endl;
	outGnu << "EOF" << endl;
    
	outGnu.close();
    
	//execute the gnu script
	string command = "chmod a+x " + gnu + "; " + gnu;
    //LogManager::log << "    Creating jpeg and eps plots in " << photometryPlotsDir << "/" << prefix  ;
    //GlobalVariables::logManager.LogManagerShowLog();    
	//system(command.c_str());
    
}
//==============================================================================






//==============================================================================
/**
 * This function writes the sub-pixels PSF mask to an ASCII file.
 * @param fileName File name
 */
void ProcessingPhotometry::processingPhotometryWritePsfToASCIIFile(string fileName)
{
    
    	ofstream out(fileName.c_str());
    	if (!out.is_open())
    	{
    		cerr << "\nError (ProcessingPhotometry::processingPhotometryWritePsfToASCIIFile): Unable to open output PSF file";
    		exit(1);
    	}
    	out.precision(10);
    
    
    	for (int i = 0; i < psfMap.rows(); i++)
    	{
    		for (int j = 0; j < psfMap.cols(); j++)
            {
    			out << (i - (psfMap.rows() - 1.) / 2.) / double(subPixelsPerPixel) << " " << (j - (psfMap.cols() - 1.) / 2.)
                       / double(subPixelsPerPixel) << " " << psfMap(i, j) << " " << endl;
    		}
            out << endl;
    	}
}
//==============================================================================




//==============================================================================
/**
 * This function writes the photometry processing info to an ASCII file.
 * @param fileName File name to be written.
 */
void ProcessingPhotometry::processingPhotometryWriteInfo(string fileName)
{
	ofstream out(fileName.c_str());
	if (!out.is_open())
	{
		cerr << "\nError (ProcessingPhotometry::processingPhotometryWriteInfo()): Unable to open output information file " << fileName;
		exit(1);
	}
    
	string baseName = outputPath + "/" + prefix;
	baseName.erase(baseName.end() - 1);
    
	out << "-----------------------------------------------------------------------------" << endl;
	out << "Information about the photometry files in this directory" << endl;
	if (photometryMethod == "WM")   //If photometry method is Weighted Mask
		out << "PHOTOMETRY METHOD: WEIGHTED MASKS" << endl;
	else if (photometryMethod == "AP") //If photometry method is Aperture
		out << "PHOTOMETRY METHOD: APERTURE" << endl;
    
	if (flatfieldCorrection)
		out << "CORRECTED FOR FLATFIELD" << endl;
	else
		out << "NOT CORRECTED FOR FLATFIELD" << endl;
    
	if (frameTransferSmearingCorrection)
		out << "CORRECTED FOR TRAILING" << endl;
	else
		out << "NOT CORRECTED FOR TRAILING" << endl;
    
	out << "Number of Telescopes (N): " << photometryNumTelescopes << endl;
	out << "-----------------------------------------------------------------------------" << endl;
	out << baseName + ".info" << "  ............. Input parameters for simulation" << endl;
	out << baseName + "*.fits" << "  ............. Raw sub-images" << endl;
	out << baseName + "_flat.fits" << "  ............. Flatfield (illumination flux 30000 e)" << endl;
	out << baseName + "_reduced0.fits" << "  ............. First image corrected by trailing, flatfield, gain and exposure (flux/s)"
    << endl;
	out << baseName + "_masks0.fits" << " ................ Masks of the first image" << endl;
	out << "Normalized flux is computed from the theoretical flux of the model PSF at the sub-pixel position of the star" << endl;
	/*
	 out << baseName + "*.phot" << " ............. Photometry of the stars on each image" << endl;
	 out << "                                      Column description:" << endl;
	 out << "                                      1: X-coordinate of star on sub-image" << endl;
	 out << "                                      2: Y-coordinate of star on sub-image" << endl;
	 out << "                                      3: Input magnitude" << endl;
	 out << "                                      4: Input flux" << endl;
	 out << "                                      5: Measured magnitude" << endl;
	 out << "                                      6: Measured normalized magnitude" << endl;
	 out << "                                      7: Measured flux" << endl;
	 out << "                                      8: Measured normalized flux" << endl;
	 out << "                                      9: Background flux" << endl;
	 */
	out << "star*.phot" << " ............. Lightcurve of each star" << endl;
	out << "                                      Column description:" << endl;
	out << "                                      1: time [s]" << endl;
	out << "                                      2: Input magnitude" << endl;
	out << "                                      3: Input flux" << endl;
	out << "                                      4: Measured magnitude" << endl;
	out << "                                      5: Measured normalized magnitude" << endl;
	out << "                                      6: Measured flux" << endl;
	out << "                                      7: Measured normalized flux" << endl;
	out << "                                      8: Background flux" << endl;
    
	out << "star_coordinates.dat" << " ............. Coordinates and magnitude of measured stars" << endl;
	out << "                                      Column description:" << endl;
	out << "                                      1: File name with lightcurve" << endl;
	out << "                                      2: X-coordinate of star on sub-image" << endl;
	out << "                                      3: Y-coordinate of star on sub-image" << endl;
	out << "                                      4: Input magnitude" << endl;
	out << "                                      5: Input flux" << endl;
    
	out << baseName + ".stat" << " ............. Statistics for all lightcurves" << endl;
	out << "                                      Column description:" << endl;
	out << "                                      1: Star file name" << endl;
	out << "                                      2: Input magnitude" << endl;
	out << "                                      3: Input flux" << endl;
	out << "                                      4: Mean input magnitude (Monte-Carlo)" << endl;
	out << "                                      5: Std. dev of input magnitude (Monte-Carlo)" << endl;
	out << "                                      6: Mean measured magnitude" << endl;
	out << "                                      7: Mean measured normalized magnitude" << endl;
	out << "                                      8: Std. dev of measured magnitude" << endl;
	out << "                                      9: RMS of measured magnitude" << endl;
	out << "                                      10: Mean measured flux" << endl;
	out << "                                      11: Mean measured normalized flux" << endl;
	out << "                                      12: Std. dev of measured flux " << endl;
	out << "                                      13: Theoretical photon noise (flux)" << endl;
	out << "                                      14: Measured noise ppm/hour for 1 telescope" << endl;
	out << "                                      15: Theoretical photon noise ppm/hour for 1 telescope" << endl;
	out << "                                      16: Measured noise ppm/hour for N telescopes" << endl;
	out << "                                      17: Theoretical photon noise ppm/hour for N telescopes" << endl;
	out << "                                      18: Measured noise ppm/hour for 8 telescopes" << endl;
	out << "                                      19: Theoretical photon noise ppm/hour for 8 telescopes" << endl;
	out << "                                      20: Measured noise ppm/hour for 16 telescopes" << endl;
	out << "                                      21: Theoretical photon noise ppm/hour for 16 telescopes" << endl;
	out << "                                      22: Measured noise ppm/hour for 24 telescopes" << endl;
	out << "                                      23: Theoretical photon noise ppm/hour for 24 telescopes" << endl;
	out << "                                      24: Measured noise ppm/hour for "<<photometryNumTelescopes<<" telescopes" << endl;
	out << "                                      25: Theoretical photon noise ppm/hour for "<<photometryNumTelescopes<<" telescopes" << endl;
    
	out << baseName + "_psf.dat" << " ............. Input PSF of image" << endl;
	out << "                                      Column description:" << endl;
	out << "                                      1: X-coordinate" << endl;
	out << "                                      2: Y-coordinate" << endl;
	out << "                                      3: Normalized flux" << endl;
    
	out << "crowding.dat" << " ............. Crowding statistics" << endl;
	out << "                                      Column description:" << endl;
	out << "                                      1: Center of range (mag)" << endl;
	out << "                                      2: Width of range (mag)" << endl;
	out << "                                      3: Total number of stars" << endl;
	out << "                                      4: Number of stars in range" << endl;
	out << "                                      5: Number of stars in range with measured mag 0.113 smaller than input mag" << endl;
	out << "                                      6: Percentage of stars in range with measured mag 0.113 smaller than input mag" << endl;
	out.close();
}
//==============================================================================




