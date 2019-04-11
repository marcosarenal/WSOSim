///////////////////////////////////////////////////////////
//  StepConvolvePSF.cpp
//  Implementation of the Class StepConvolvePSF
//  Created on:      23-Oct-2012 2:00:02 PM
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



#include "StepConvolvePSF.h"



//==============================================================================
//==============================================================================
/**
 * Constructor method
 */
StepConvolvePSF::StepConvolvePSF(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
StepConvolvePSF::~StepConvolvePSF(){}
//==============================================================================



//==============================================================================
//Functions:
//==============================================================================
/**
 * This function convolves the PSF.
 */
void StepConvolvePSF::StepConvolvePSFapplication(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    convolutionMethod = p_DataSet->datasetGetconvolutionMethod();
    threads = p_DataSet->datasetGetthreads();
    useFFTWisdom = p_DataSet->datasetGetuseFFTWisdom();
    fftWisdomPath = p_DataSet->datasetGetfftWisdomPath();
    memoryLimit = p_DataSet->datasetGetmemoryLimit();
    
    LogManager::log << "    Convolving with PSF...";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //Retrieving the PSF Map from DataSet
//    m_DataSet.datasetGetPSFMap(psfMap);

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

    
    //completing the convolved map with the subpixel map
    convolvedMap.resize(subPixelMap.shape());
    convolvedMap = subPixelMap;
    
    
    //Checking what type of convolution is required
    if (convolutionMethod == "FFT")
    {
        
        StepConvolvePSF::StepConvolvePSFcomputeFFT();
    }
    else
    {
        StepConvolvePSF::StepConvolvePSFcomputeReal(); 
    }
    if (max(convolvedMap) == 0)
    {

        std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFapplication()): The simulated image does not contain any flux." << std::endl;
        exit(1);
    }
        
    //Setting the new calculated pixelMap into the DataSet
    m_DataSet.datasetSetpixelMap(pixelMap);
    
    //Set the calculated convolved map as the subPixelMap    
    m_DataSet.datasetSetsubPixelMap(convolvedMap);

    //Free memory
    convolvedMap.free();
}
//==============================================================================




//==============================================================================
/**
 * This function computes the convolution process to the psfMap array.
 */
void StepConvolvePSF::StepConvolvePSFcomputeFFT()
{
    LogManager::log <<"    FFT convolution";
    GlobalVariables::logManager.LogManagerShowLog();
    
    //Check whether PSF is empty
    if (psfMap.rows() == 0)
    {
        std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeFFT()): Mask for FFT convolution is empty." << std::endl;
        exit(1);
    }
     
    //Check if psfMap is quadratic
    if (psfMap.rows() != psfMap.cols())
	{
		std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeFFT()): PSF must be quadratic." << std::endl;
        exit(1);
	}

    //Check if psfMap has odd size
    if (psfMap.rows() % 2 == 0)
    {
        std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeFFT()): PSF must have odd size." << std::endl;
        exit(1);
    }
    
    
    int psfSize = floor(psfMap.rows() / 2);
    
    //Calculation of the amount of memory required for the subPixelMap, in bytes
    float subPixelMapMemory;
    subPixelMapMemory = 4*(subPixelMap.rows() + psfSize) * 4*(subPixelMap.cols() + psfSize); //int = 4 bytes, float = 4 bytes
    
    //If there is enough memory available
    if (subPixelMapMemory < memoryLimit)
    {
        StepConvolvePSF::StepConvolvePSFconvolveFFT(subPixelMap, psfMap, convolvedMap, threads, useFFTWisdom, fftWisdomPath); //if image is too big split it up to avoid memory overflow
    }
    else
    {
        LogManager::log <<"    Splitting array for convolution";
        GlobalVariables::logManager.LogManagerShowLog();

        Array<float, 2> tempConv;
        tempConv = 0.0;

        int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        int i = 0, j = 0;
        int upperColsLimit = int(sqrt(double(memoryLimit)) - double(psfMap.rows()) / 2.);
        int n = 1, maxCols;

        //set size of columns and rows of the segment that is FFTed such that the FFT receives a field 2^n by 2^n (much faster FFT!)
        //at image edges this fails
    while (pow(2.0, n) - 3 * psfSize <= upperColsLimit)
    {
        n++;
    }

            maxCols = pow(2.0, n - 1) - 3 * psfSize - 1;

        
        LogManager::log <<"    Convoluting ...";
        GlobalVariables::logManager.LogManagerShowLog();
        
		int counter = 0;
		double maxCounter = 4. + 1. * (subPixelMap.rows() - psfSize) * (subPixelMap.cols() - psfSize) / pow2(1. * maxCols);
		while ((j) * maxCols + psfSize < subPixelMap.rows())
		{
			i = 0;
			y0 = j * maxCols;
			y1 = (j + 1) * maxCols + 2 * psfSize;
            
			if (y1 >= subPixelMap.rows())
            {
				y1 = subPixelMap.rows() - 1;
            }
            
            
			while ((i) * maxCols + psfSize < subPixelMap.cols())
			{
                //LogManager::log <<"    "<< floor((counter++) * 100. / maxCounter )<< " % ... ";
                //GlobalVariables::logManager.LogManagerShowLog();
				x0 = i * maxCols;
				x1 = (i + 1) * maxCols + 2 * psfSize;
                
				if (x1 >= subPixelMap.cols())
                {
					x1 = subPixelMap.cols() - 1;
                }
				// if ( x0 == 0 ) convolvedMap ( Range ( x0 , x1 - psfSize ), Range::all() ) = tempConv ( Range ( 0, toEnd - psfSize ), Range::all() );
				//      else if ( x1 == subPixelMap.cols() - 1 ) convolvedMap ( Range ( x0 + psfSize, x1 ), Range::all() ) = tempConv ( Range ( psfSize, toEnd ), Range::all() );
				//    else convolvedMap ( Range ( x0 + psfSize, x1 - psfSize ), Range::all() ) = tempConv ( Range ( psfSize, toEnd - psfSize ), Range::all() );
                
                //Call to the PSF convolution function
				StepConvolvePSF::StepConvolvePSFconvolveFFT(subPixelMap(Range(x0, x1), Range(y0, y1)), psfMap, tempConv, threads, useFFTWisdom, fftWisdomPath);
                
                
				if (x0 == 0 && y0 == 0)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0, y1 - psfSize)) = tempConv(Range(0, toEnd - psfSize), Range(0, toEnd - psfSize));
                }
				else if (x0 == 0 && y0 > 0 && y1 < subPixelMap.rows() - 1)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(0, toEnd - psfSize), Range(psfSize, toEnd - psfSize));
                }
				else if (x0 == 0 && y1 == subPixelMap.rows() - 1)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0 + psfSize, y1)) = tempConv(Range(0, toEnd - psfSize), Range(psfSize, toEnd));   
                }
				else if (x0 > 0 && y0 == 0 && x1 < subPixelMap.cols() - 1)
                {
                    convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0, y1 - psfSize)) = tempConv(Range(psfSize, toEnd - psfSize), Range(0, toEnd - psfSize));
                }
				else if (x1 == subPixelMap.cols() - 1 && y0 == 0)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0, y1 - psfSize)) = tempConv(Range(psfSize, toEnd), Range(0, toEnd - psfSize));
                }
				else if (x1 == subPixelMap.cols() - 1 && y1 < subPixelMap.rows() - 1)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(psfSize, toEnd), Range(psfSize, toEnd - psfSize));
                }
				else if (x1 < subPixelMap.cols() - 1 && y1 == subPixelMap.rows() - 1)
                {
					convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0 + psfSize, y1)) = tempConv(Range(psfSize, toEnd - psfSize), Range(psfSize, toEnd));
				}
                else if (x1 == subPixelMap.cols() - 1 && y1 == subPixelMap.rows() - 1)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0 + psfSize, y1)) = tempConv(Range(psfSize, toEnd), Range(psfSize, toEnd));
				}
                else
				{
                    convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(psfSize, toEnd - psfSize), Range(psfSize, toEnd - psfSize));
				}
                i++;
			//	tempConv.free();
			}
			j++;
		}
            
		//Free temporary conv. map		
        tempConv.free();
		
	}
    
}
//==============================================================================




//==============================================================================
/**
 * This function computes the real convolution process to the psfMap array.
 */
void StepConvolvePSF::StepConvolvePSFcomputeReal()
{
	LogManager::log <<"    Real convolution";
    GlobalVariables::logManager.LogManagerShowLog();
    
	//Check whether PSF is empty
    if (psfMap.rows() == 0)
	{
		std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeReal()): Mask for FFT convolution is empty." << std::endl;
        exit(1);
	}
       
	//Check if psfMap is quadratic
    if (psfMap.rows() != psfMap.cols())
	{
		std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeReal()): PSF must be quadratic." << std::endl;
        exit(1);
	}

    //Check if psfMap has odd size
	if (psfMap.rows() % 2 == 0)
	{
		std::cerr << "\nError (StepConvolvePSF::StepConvolvePSFcomputeReal()): PSF must have odd size." << std::endl;
        exit(1);
	}
     
    
  	int psfSize = floor(psfMap.rows() / 2);
    
    //Calculation of the amount of memory required for the subPixelMap, in bytes
    float subPixelMapMemory;
    subPixelMapMemory = 4*(subPixelMap.rows() + psfSize) * 4*(subPixelMap.cols() + psfSize); //int = 4 bytes, float = 4 bytes
    
    //If there is enough memory available
	if (subPixelMapMemory < memoryLimit)
    {
		StepConvolvePSF::StepConvolvePSFconvolveReal(subPixelMap, psfMap, convolvedMap); //if image is too big split it up to avoid memory overflow
	}
    else
	{
		LogManager::log <<"    Splitting array for convolution";
        GlobalVariables::logManager.LogManagerShowLog();
        
		Array<float, 2> tempConv;
		int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
		int i = 0, j = 0;
		int upperColsLimit = int(sqrt(double(memoryLimit)) - double(psfMap.rows()) / 2.);
		int n = 1, maxCols;
        
		//set size of columns and rows of the segment that is FFTed such that the FFT receives a field 2^n by 2^n (much faster FFT!)
		//at image edges this fails
		while (pow(2.0, n) - 3 * psfSize <= upperColsLimit)
        {
			n++;
        }
        
		maxCols = pow(2.0, n - 1) - 3 * psfSize - 1;
        
        LogManager::log <<"    Convoluting ...";
        GlobalVariables::logManager.LogManagerShowLog();
        
		int counter = 0;
		double maxCounter = 4. + 1. * (subPixelMap.rows() - psfSize) * (subPixelMap.cols() - psfSize) / pow2(1. * maxCols);
		while ((j) * maxCols + psfSize < subPixelMap.rows())
		{
			i = 0;
			y0 = j * maxCols;
			y1 = (j + 1) * maxCols + 2 * psfSize;
            
			if (y1 >= subPixelMap.rows())
            {
				y1 = subPixelMap.rows() - 1;
            }
            
			while ((i) * maxCols + psfSize < subPixelMap.cols())
			{
                //LogManager::log <<"    "<< floor((counter++) * 100. / maxCounter )<< " % ... ";
                //GlobalVariables::logManager.LogManagerShowLog();
				x0 = i * maxCols;
				x1 = (i + 1) * maxCols + 2 * psfSize;
                
				if (x1 >= subPixelMap.cols())
					x1 = subPixelMap.cols() - 1;
				// if ( x0 == 0 ) convolvedMap ( Range ( x0 , x1 - psfSize ), Range::all() ) = tempConv ( Range ( 0, toEnd - psfSize ), Range::all() );
				//      else if ( x1 == subPixelMap.cols() - 1 ) convolvedMap ( Range ( x0 + psfSize, x1 ), Range::all() ) = tempConv ( Range ( psfSize, toEnd ), Range::all() );
				//    else convolvedMap ( Range ( x0 + psfSize, x1 - psfSize ), Range::all() ) = tempConv ( Range ( psfSize, toEnd - psfSize ), Range::all() );
                
				StepConvolvePSF::StepConvolvePSFconvolveReal(subPixelMap(Range(x0, x1), Range(y0, y1)), psfMap, tempConv);
                
                
                
				if (x0 == 0 && y0 == 0)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0, y1 - psfSize)) = tempConv(Range(0, toEnd - psfSize), Range(0, toEnd - psfSize));
                }
				else if (x0 == 0 && y0 > 0 && y1 < subPixelMap.rows() - 1)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(0, toEnd - psfSize), Range(psfSize, toEnd - psfSize));
                }
				else if (x0 == 0 && y1 == subPixelMap.rows() - 1)
                {
                    convolvedMap(Range(x0, x1 - psfSize), Range(y0 + psfSize, y1)) = tempConv(Range(0, toEnd - psfSize), Range(psfSize, toEnd));   
                }
				else if (x0 > 0 && y0 == 0 && x1 < subPixelMap.cols() - 1)
                {
                    convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0, y1 - psfSize)) = tempConv(Range(psfSize, toEnd - psfSize), Range(0, toEnd - psfSize));
                }
				else if (x1 == subPixelMap.cols() - 1 && y0 == 0)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0, y1 - psfSize)) = tempConv(Range(psfSize, toEnd), Range(0, toEnd - psfSize));
                }
				else if (x1 == subPixelMap.cols() - 1 && y1 < subPixelMap.rows() - 1)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(psfSize, toEnd), Range(psfSize, toEnd - psfSize));
                }
				else if (x1 < subPixelMap.cols() - 1 && y1 == subPixelMap.rows() - 1)
                {
					convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0 + psfSize, y1)) = tempConv(Range(psfSize, toEnd - psfSize), Range(psfSize, toEnd));
				}
                else if (x1 == subPixelMap.cols() - 1 && y1 == subPixelMap.rows() - 1)
				{
                    convolvedMap(Range(x0 + psfSize, x1), Range(y0 + psfSize, y1)) = tempConv(Range(psfSize, toEnd), Range(psfSize, toEnd));
				}
                else
				{
                    convolvedMap(Range(x0 + psfSize, x1 - psfSize), Range(y0 + psfSize, y1 - psfSize)) = tempConv(Range(psfSize, toEnd - psfSize), Range(psfSize, toEnd - psfSize));
				}
				i++;
				tempConv.free();
			}
			j++;
		}
	}
    
}
//==============================================================================




//==============================================================================
/**
 *
 * @param tempSubPixelMap
 * @param tempPsfMap
 * @param convolvedMap
 * @param threads
 * @param useFFTWisdom
 * @param fftWisdomPath
 */
void StepConvolvePSF::StepConvolvePSFconvolveFFT(Array<float, 2> tempSubPixelMap, Array<float, 2> tempPsfMap, Array<float, 2> &convolvedMap, int threads, bool useFFTWisdom,
                                              std::string fftWisdomPath)
{
    
	int numRows = tempSubPixelMap.rows();
	int numCols = tempSubPixelMap.cols();
    
	StepConvolvePSF::StepConvolvePsfFftRemapArrays(tempSubPixelMap, tempPsfMap);
	int numColsRemapped = tempSubPixelMap.cols();
    
	Array<complex<float> , 2> fdata, fmask;
    
	//compute FFT of tempSubPixelMap and tempPsfMap
	FFT::forward(tempSubPixelMap, fdata, threads, useFFTWisdom, fftWisdomPath);
	tempSubPixelMap.free();
    
	FFT::forward(tempPsfMap, fmask, threads, useFFTWisdom, fftWisdomPath);
	tempPsfMap.free();
    
	//multiply the transformed
	Array<complex<float> , 2> product(fdata.shape());
	product = fdata * fmask;
    
	fdata.free();
	fmask.free();
    
	//transform product back to real space
	FFT::backward(numColsRemapped, product, convolvedMap, threads, useFFTWisdom, fftWisdomPath);
	
    //normalize
	convolvedMap /= (convolvedMap.rows() * convolvedMap.cols());
    
	convolvedMap.resizeAndPreserve(numRows, numCols);
    
	product.free();
    
}
//==============================================================================




//==============================================================================
/**
 *
 * @param tempSubPixelMap
 * @param tempPsfMap
 * @param convolvedMap
 */
void StepConvolvePSF::StepConvolvePSFconvolveReal(Array<float, 2> tempSubPixelMap, Array<float, 2> tempPsfMap, Array<float, 2> &convolvedMap)
{
	if (tempPsfMap.rows() == 0)
	{
		convolvedMap.resize(tempSubPixelMap.shape());
		convolvedMap = tempSubPixelMap;
		return;
	}
    
	if (tempPsfMap.rows() != tempPsfMap.cols())
		return; //tempPsfMap must be quadratic
	if (tempPsfMap.rows() % 2 == 0)
		return; //tempPsfMap must have odd size
    
	int numRows = tempSubPixelMap.rows();
	int numCols = tempSubPixelMap.cols();
	int maskCenter = (tempPsfMap.rows() - 1) / 2;
    
	StepConvolvePSF::StepConvolvePsfRealRemapArrays(tempSubPixelMap, tempPsfMap);
    
	convolvedMap.resize(tempSubPixelMap.shape());
	convolvedMap = 0.0;
    
	for (int i = maskCenter; i < tempSubPixelMap.rows() - maskCenter; i++)
		for (int j = maskCenter; j < tempSubPixelMap.cols() - maskCenter; j++)
			if (tempSubPixelMap(i, j) > 0.)
			{
				convolvedMap(Range(i - maskCenter, i + maskCenter), Range(j - maskCenter, j + maskCenter)) += (tempSubPixelMap(i, j) * tempPsfMap);
			}
    
	convolvedMap(Range(0, numRows - 1), Range(0, numCols - 1)) = convolvedMap(Range(maskCenter, convolvedMap.rows() - 1 - maskCenter), Range(maskCenter,
                                                                                                                                             convolvedMap.cols() - 1 - maskCenter));
    
	convolvedMap.resizeAndPreserve(numRows, numCols);
    
}
//==============================================================================




//==============================================================================
/**
 * This function remaps the subPixelMap and psfMap for FFT convolution method
 * The subPixelReMap array is extended to the size (subPixelMap.size+psfMap.size) to avoid problems when convoluting
 *
 * @param subPixelReMap
 * @param psfReMap
 */
void StepConvolvePSF::StepConvolvePsfFftRemapArrays(Array<float, 2> &subPixelReMap, Array<float, 2> &psfReMap)
{
	//remap the subPixelMap and psfMap
	//the subPixelReMap array is extended to the size (subPixelMap.size+psfMap.size) to avoid problems when convoluting
	//the psfReMap array is extended to the same size. the values of the psfReMap are copied into the corners of this array (newMask)
	int drows0 = subPixelReMap.rows();
	int dcols0 = subPixelReMap.cols();
    
	subPixelReMap.resizeAndPreserve(drows0 + floor(psfReMap.rows() / 2), dcols0 + floor(psfReMap.cols() / 2));
	subPixelReMap(Range(drows0, toEnd), Range::all()) = 0.;
	subPixelReMap(Range::all(), Range(dcols0, toEnd)) = 0.;
    
	Array<float, 2> newMask(subPixelReMap.shape());
	newMask = 0.;
	int middle = (psfReMap.rows() - 1) / 2.; //index of middle pixel of psfReMap
    
	//copy lower right part of psfReMap into upper left corner of newMask
	newMask(Range(fromStart, middle), Range(fromStart, middle)) = psfReMap(Range(middle, toEnd), Range(middle, toEnd));
	//copy upper right part of psfReMap into lower left corner of newMask
	newMask(Range(newMask.rows() - middle, toEnd), Range(fromStart, middle)) = psfReMap(Range(fromStart, middle - 1), Range(middle, toEnd));
	//copy lower left part of psfReMap into upper right corner of newMask
	newMask(Range(fromStart, middle), Range(newMask.cols() - middle, toEnd)) = psfReMap(Range(middle, toEnd), Range(fromStart, middle - 1));
	//copy upper left part of psfReMap into lower right corner of newMask
	newMask(Range(newMask.rows() - middle, toEnd), Range(newMask.cols() - middle, toEnd)) = psfReMap(Range(fromStart, middle - 1), Range(
                                                                                                                                       fromStart, middle - 1));
    
	psfReMap.resize(newMask.shape());
	psfReMap = newMask;
    
	newMask.free();
}
//==============================================================================




//==============================================================================
/**
 * This function remaps the subPixelMap and psfMap for Real convolution method
 * The subPixelReMap array is extended to the size (subPixelMap.size+psfMap.size) to avoid problems when convoluting
 * 
 * @param subPixelReMap
 * @param psfReMap
 */
void StepConvolvePSF::StepConvolvePsfRealRemapArrays(Array<float, 2> &subPixelReMap, Array<float, 2> &psfReMap)
{
	//remap the subPixelMap, psfMap is unmodified
	//the subPixelReMap array is extended to the size (subPixelMap.size+psfMap.size) to avoid problems when convolving
	int drows0 = subPixelReMap.rows();
	int dcols0 = subPixelReMap.cols();
    
	subPixelReMap.resizeAndPreserve(drows0 + psfReMap.rows(), dcols0 + psfReMap.cols());
    
	Array<float, 2> newData(subPixelReMap.shape());
    
	newData = 0.;
	int halfMaskRows = floor(psfReMap.rows() / 2);
	int halfMaskCols = floor(psfReMap.cols() / 2);
    
	newData(Range(halfMaskRows, newData.rows() - halfMaskRows), Range(halfMaskCols, newData.cols() - halfMaskCols)) = subPixelReMap(Range(fromStart,
                                                                                                                                        subPixelReMap.rows() - psfReMap.rows() - 1), 
                                                                                                                                    Range(fromStart, 
                                                                                                                                        subPixelReMap.cols() - psfReMap.cols() - 1));
    
	subPixelReMap = newData;
	newData.free();
}
//==============================================================================

