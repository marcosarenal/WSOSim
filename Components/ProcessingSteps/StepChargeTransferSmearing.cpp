///////////////////////////////////////////////////////////
//  StepChargeTransferSmearing.cpp
//  Implementation of the Class StepChargeTransferSmearing
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



#include "StepChargeTransferSmearing.h"
#include "MathTools.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepChargeTransferSmearing::StepChargeTransferSmearing(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepChargeTransferSmearing::~StepChargeTransferSmearing(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function adds the Charge Transfer Smearing effect to the image by
 * summing up the intensity of all stars in the same columns as the sub image. 
 * It must be taken into account that the readout direction is y.
 * This function also generates an smearingMap; this is the image to be obtained in 
 * the over-scan. This smearingMap (over-scan) is used in the Photometry to remove 
 * the smearing effect from the real image (pixelMap).
 */
void StepChargeTransferSmearing::StepChargeTransferSmearingapplication(DataSet &m_DataSet)
{

        //Pointing to the DataSet
        p_DataSet = &m_DataSet;
        	
        
        //Retrieving parameters from DataSet
        exposureTime = p_DataSet->datasetGetexposureTime();
        readOutTime = p_DataSet->datasetGetreadOutTime();
        numSmearingOverscanRows = p_DataSet->datasetGetnumSmearingOverscanRows();
        subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
        subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();        
        subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
        background = p_DataSet->datasetGetbackground();

        subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
        
        //Retrieving the PSF Map from DataSet
//        m_DataSet.datasetGetPSFMap(psfMap);
        
        //Initialize PSF map  
        psfMap.resize(m_DataSet.datasetGetPSFMap().extent(0),m_DataSet.datasetGetPSFMap().extent(1));
        psfMap = 0.0;
        
        //Retrieving the PSF Map from DataSet
        psfMap = m_DataSet.datasetGetPSFMap();
        
        //Initializing and retrieving the pixel map from DataSet
        pixelMap.resize(m_DataSet.datasetGetpixelMap().extent(0),m_DataSet.datasetGetpixelMap().extent(1));
        pixelMap = 0.0;
        pixelMap = m_DataSet.datasetGetpixelMap();
        
        //Initializing and retrieving the subpixel map from DataSet
        subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
        subPixelMap = 0.0;
        subPixelMap = m_DataSet.datasetGetsubPixelMap();
        
                 
        //Retrieving the starListOnCCD contains for each star, its X and Y position in CCD, magnitude, RA, 
        //declination and identification number.
        starListOnCCD.resize(m_DataSet.datasetGetStarListOnCCD().extent(0),m_DataSet.datasetGetStarListOnCCD().extent(1));
        starListOnCCD = 0.0;
        starListOnCCD = m_DataSet.datasetGetStarListOnCCD();
  
        //Temporary Blitz array to generate a smearing map              
        Array<float, 1> conv;
        conv = 0.0;


        //Sum up the intensity of all stars in the same columns as the sub image.
        StepChargeTransferSmearing::StepChargeTransferSmearinggenerateSmearingMap(conv);
               
          

        //spread the flux across all rows at subPixelMap
        for (unsigned int i = 0; i < conv.size(); i++)
        {
            //The smearing is added to the subPixelMap
            subPixelMap(i, Range::all()) += conv(i);
        }

                
        //Compute the trailing at pixel level (for photometry)
        smearingMap.resize(subFieldSizeY, numSmearingOverscanRows);
        smearingMap = 0.0;                
                
        //Rebinnig from subpixel to pixel level and added smearing values to smearingMap
        for (int i = 0; i < subFieldSizeX; i++)
        {
            smearingMap(i, Range::all()) = subPixelsPerPixel * (sum(conv(Range(i * subPixelsPerPixel, (i + 1) * subPixelsPerPixel - 1))));
        }
         
        //Add the background level. Since the background is spread across complete CCD,
        //its contribution is simply transfer time * background
        smearingMap += background * readOutTime;

        //Free the blitz array
        conv.free();

        //Set the empty pixel map and subpixel map.
        m_DataSet.datasetSetsmearingMap(smearingMap);
        m_DataSet.datasetSetsubPixelMap(subPixelMap);

    

        LogManager::log <<"    Successfully added Charge Transfer Smearing effect.";
        GlobalVariables::logManager.LogManagerShowLog();     
	

}
//==============================================================================




//==============================================================================
/**
 * Sum up the intensity of all stars in the same columns as the sub image. 
 * It must be taken into account that the readout direction is y. 
 * @return &conv is the smearing map generated to be added to the subPixel map.
 */
void StepChargeTransferSmearing::StepChargeTransferSmearinggenerateSmearingMap(Array<float, 1> &conv)
{

    //Retrieve parameters from Dataset
    fluxm0 = p_DataSet->datasetGetfluxm0();
    areaTelescope = p_DataSet->datasetGetareaTelescope();
    transEff = p_DataSet->datasetGettransEff();
    quantEff = p_DataSet->datasetGetquantEff();
    ccdSizeY = p_DataSet->datasetGetccdSizeY();
    edgePixels = p_DataSet->datasetGetedgePixels();
    
    //Defined the factor of conversion from magnitudes to flux
    double mag2fluxFactor = fluxm0 * exposureTime * areaTelescope * transEff * quantEff;
    
    //Define a blitz array to sum the flux in one column (y-direction) at subpixel level
    Array<float, 1> sumFlux(subFieldSizeX*subPixelsPerPixel);
    sumFlux = 0.0;
    
    //Iterator on the x coordinate
    int xCoord;
    
          
    //Take every source in the starlistOnCCD array and sum the fluxes of each row
    for (int detectedStar = 0; detectedStar < starListOnCCD.rows(); detectedStar++)
    {
        //compute the x coordinate of the stars at intra pixel level
        xCoord = int(round(starListOnCCD(detectedStar, 0) - subFieldZeroX * subPixelsPerPixel));

        if (xCoord >= 0 && xCoord < subFieldSizeX*subPixelsPerPixel)
        {
            sumFlux(xCoord) += mag2fluxFactor * pow(10., -0.4 * starListOnCCD(detectedStar, 2));

        }
    }


    //Do the same for the PSF
    Array<float, 1> sumPSF(psfMap.cols());
    for (int detectedStar = 0; (unsigned)detectedStar < sumPSF.size(); detectedStar++)
    {
        sumPSF(detectedStar) = sum(psfMap(detectedStar, Range::all()));
    }
        
    //convolve the two arrays
    conv.resize(sumFlux.rows());
    conv = 0.0;
    StepChargeTransferSmearing::StepChargeTransferSmearingconvolve(sumFlux, sumPSF, conv);

    //scale conv with the number of pixels and the transfer time
    double factor = (readOutTime / (exposureTime * (ccdSizeY - edgePixels * 2) * subPixelsPerPixel));
    conv *= factor;

    sumFlux.free();
    sumPSF.free();
}
//==============================================================================



//==============================================================================
/**
 * This function takes two 1-D arrays and returns the convolved 1-D array termed conv.
 * @params data First input array. 
 * @params mask Second input array. 
 */
void StepChargeTransferSmearing::StepChargeTransferSmearingconvolve(Array<float, 1> data, Array<float, 1> mask, Array<float, 1> &conv)
{

    if (mask.size() == 0)
    {
        conv.resize(data.shape());
        conv = data;
        return;
    }

    if (mask.size() % 2 == 0)
    return; //mask must have odd size

    int size = data.size();

    StepChargeTransferSmearing::StepChargeTransferSmearingremapArrays(data, mask);

    Array<complex<float> , 1> fdata, fmask, product(data.shape());

    //compute FFT of data and mask
    FFT::forward(data, fdata);
    FFT::forward(mask, fmask);

    //multiply the transformed
    product = fdata * fmask;

    //transform product back to real space
    FFT::backward(product, conv);

    //normalize
    conv /= (conv.size()); 

    conv.resizeAndPreserve(size);

    fdata.free();
    fmask.free();
    product.free();
}
//==============================================================================



//==============================================================================
/**
 * This function extends both input arrays to the same size in order to avoid later 
 * problems when performing convolution.
 * @params data First input array. 
 * @params mask Second input array. 
 */
void StepChargeTransferSmearing::StepChargeTransferSmearingremapArrays(Array<float, 1> &data, Array<float, 1> &mask)
{   
    //Remap the data and mask
    //The data array is extended to the size (data.size+mask.size) to avoid problems when convolving
    //The mask array is extended to the same size. The values of the mask are copied into the corners of this array (newMask)
    int size0 = data.size();

    data.resizeAndPreserve(size0 + floor(mask.size() / 2));
    data(Range(size0, toEnd)) = 0.;

    Array<float, 1> newMask(data.shape());
    newMask = 0.;
    int middle = (mask.size() - 1) / 2.; //index of middle pixel of mask

    //copy right part of mask into left corner of newMask
    newMask(Range(fromStart, middle)) = mask(Range(middle, toEnd));
    //copy left part of mask into right corner of newMask
    newMask(Range(newMask.size() - middle, toEnd)) = mask(Range(fromStart, middle - 1));

    mask.resize(newMask.shape());
    mask = newMask;

    newMask.free();
}
//==============================================================================
