///////////////////////////////////////////////////////////
//  StepCosmics.cpp
//  Implementation of the Class StepCosmics
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



#include "StepCosmics.h"





//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepCosmics::StepCosmics(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepCosmics::~StepCosmics(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 *This function applies cosmic hit to the image.
 */
void StepCosmics::StepCosmicsapplication(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;

    //Retrieving parameters from DataSet
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    fullWellSat = p_DataSet->datasetGetfullWellSat();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    exposureTime = p_DataSet->datasetGetexposureTime();
    pixelSize = p_DataSet->datasetGetpixelSize();
    
    cosmicHitRate = p_DataSet->datasetGetcosmicHitRate();
    cosmicsLength = p_DataSet->datasetGetcosmicsLength();
    cosmicsWidth = p_DataSet->datasetGetcosmicsWidth();
    cosmicsSatFactor = p_DataSet->datasetGetcosmicsSatFactor();

    //Initializing and retrieving the subpixel map from DataSet
    subPixelMap.resize(m_DataSet.datasetGetsubPixelMap().extent(0),m_DataSet.datasetGetsubPixelMap().extent(1));
    subPixelMap = 0.0;
    subPixelMap = m_DataSet.datasetGetsubPixelMap();
     
    // X position of cosmic on ccd
	DiscreteUniform<int, ranlib::MersenneTwister, ranlib::independentState> randX(subPixelMap.rows() - 1);
    
    // Y position of cosmic on ccd
	DiscreteUniform<int, ranlib::MersenneTwister, ranlib::independentState> randY(subPixelMap.cols() - 1); 
    
    //Angle of cosmic axis
	Uniform<double, ranlib::MersenneTwister, ranlib::independentState> theta; 
    
    //Set the intensity as from a Normal distribution with mean value =cosmicsSatFactor* fullWellSat and standard deviation = sqrt(cosmicsSatFactor*fullWellSat).
	Normal<double, ranlib::MersenneTwister, ranlib::independentState> intensity(cosmicsSatFactor*fullWellSat, sqrt(cosmicsSatFactor*fullWellSat));

    //Set the sigmaX as from a Normal distribution with mean value = cosmicsLength*subPixelsPerPixel and standard deviation = sqrt(cosmicsLength*subPixelsPerPixel).
	Normal<double, ranlib::MersenneTwister, ranlib::independentState> sigmaX(cosmicsLength*subPixelsPerPixel, sqrt(cosmicsLength*subPixelsPerPixel));

    //Set the sigmaY as from a Normal distribution with mean value = cosmicsWidth*subPixelsPerPixel and standard deviation = sqrt(cosmicsWidth*subPixelsPerPixel).
	Normal<double, ranlib::MersenneTwister, ranlib::independentState> sigmaY(cosmicsWidth*subPixelsPerPixel, sqrt(cosmicsWidth*subPixelsPerPixel));
    
    //Seed these random parameters
	randX.seed(DataSet::seedRNG);
	randY.seed(DataSet::seedRNG + 1);
	theta.seed(DataSet::seedRNG + 2);
	intensity.seed(DataSet::seedRNG + 3);
	sigmaX.seed(DataSet::seedRNG + 4);
	sigmaY.seed(DataSet::seedRNG + 5);

	//Determine how many hits are registered during the exposure
    //compute size of sub-field in cm^2 from [pixel]*[pixel]*[microns/pixel]^2*[cm^2/microns^2]
    ccdSubFieldSizeSqcm = subFieldSizeX * subFieldSizeY * pow2(pixelSize) * 1E-8;

	//hitRate is in events/cm^2/min, so we have to convert this to a rate per seconds
	meanEvents = ccdSubFieldSizeSqcm * exposureTime * cosmicHitRate / 60.; //mean number of events during the exposure time and on the complete sub-field

	int numEvents = Statistics::getPoisson(meanEvents, DataSet::seedRNG + 5);

    //For each event in numEvents, one cosmic ray is generated
	for (int i = 0; i < numEvents; i++)
    {
		StepCosmics::StepCosmicsadd2DGauss(randX.random(), randY.random(), max(0., intensity.random()), max(1.,
				sigmaX.random()), max(1., sigmaY.random()), theta.random() * Constants::Pi2);
    }
    
    
    //Setting the new calculated subPixelMap into the DataSet
    m_DataSet.datasetSetsubPixelMap(subPixelMap);  

        
    LogManager::log <<"    Affected by "<< numEvents << " cosmic hits.";
    GlobalVariables::logManager.LogManagerShowLog();    
}
//==============================================================================




//==============================================================================
/**
 * This function is thrigered once for each event in numEvents. 
 * For each call, one cosmic ray is generated and added to the subPixelMap.
 * This function generates a 2 dimensional Gaussian  distribution taking the
 * following parameters as input:
 * @param centerX X coordinate of the cosmic ray on the CCD.
 * @param centerY Y coordinate of the cosmic ray on the CCD.
 * @param intensity Intensity of the cosmic hit.
 * @param sigmax Length of the cosmic hit.
 * @param sigmay Width of the cosmic hit.
 * @param theta Angle of cosmic axis.
 */
void StepCosmics::StepCosmicsadd2DGauss(int centerX, int centerY, double intensity, double sigmax, double sigmay, double theta)
{
    
	int size = max(double(subPixelsPerPixel), 10 * max(sigmax, sigmay));////care: factor 10 is arbitrary and set to make the mask sufficiently large
	
    //mask has at least size subPixelsPerPixel to ensure we cover a complete pixel (to compute normFactor for intensity)
	//mask must have odd size
	if (size % 2 == 0)
		size++;

	int middle = (size - 1) / 2;

	double a = pow2(cos(theta)) / 2. / pow2(sigmax) + pow2(sin(theta)) / 2. / pow2(sigmay);
	double b = -sin(2. * theta) / 4. / pow2(sigmax) + sin(2. * theta) / 4. / pow2(sigmay);
	double c = pow2(sin(theta)) / 2. / pow2(sigmax) + pow2(cos(theta)) / 2. / pow2(sigmay);

	Array<float, 2> mask(size, size);
	firstIndex i;
	secondIndex j;
    
    //Generate a Gaussian-shaped mask
	mask = exp(-(a * pow2(i - middle) + 2. * b * (i - middle) * (j - middle) + c * pow2(j - middle)));

	//normalize the mask such that the flux in the central pixel will be equivalent to the intensity given
	double normFactor = sum(mask(Range(middle - subPixelsPerPixel / 2, middle + subPixelsPerPixel / 2), Range(middle - subPixelsPerPixel
			/ 2, middle + subPixelsPerPixel / 2)));
	mask *= intensity / normFactor;

	//copy cosmic into ccd subPixelMap
	int x0, y0, x1, y1;
	//coordinates of the mask in the ccd frame (have to be cut if cosmic is at the border of the ccd)
	x0 = max(0, centerX - middle);
	x1 = min(subPixelMap.rows() - 1, centerX + middle);
	y0 = max(0, centerY - middle);
	y1 = min(subPixelMap.cols() - 1, centerY + middle);

	//corner points of the trimmed mask (to match the coordinates in the ccd above)
	int maskx0, masky0, maskx1, masky1;
	if (centerX - middle < 0)
		maskx0 = middle - centerX;
	else
		maskx0 = 0;
	if (centerX + middle > subPixelMap.rows())
		maskx1 = mask.rows() - 2 + subPixelMap.rows() - middle - centerX;
	else
		maskx1 = mask.rows() - 1;
	if (centerY - middle < 0)
		masky0 = middle - centerY;
	else
		masky0 = 0;
	if (centerY + middle > subPixelMap.cols())
		masky1 = mask.cols() - 2 + subPixelMap.cols() - middle - centerY;
	else
		masky1 = mask.cols() - 1;

    //The cosmic hit is added to the subPixelMap.
	subPixelMap(Range(x0, x1), Range(y0, y1)) += mask(Range(maskx0, maskx1), Range(masky0, masky1));
    
}
//==============================================================================
