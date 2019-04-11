///////////////////////////////////////////////////////////
//  StepConvolvePSF.h
//  Implementation of the Class StepConvolvePSF
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



#ifndef StepConvolvePSF_H_
#define StepConvolvePSF_H_

#include "DataSet.h"
#include "LogManager.h"
#include "FFT.h"
#include "FileUtilities.h"

class StepConvolvePSF
{

public:
	StepConvolvePSF();
	virtual ~StepConvolvePSF();
        void StepConvolvePSFapplication(DataSet &m_DataSet);
        void StepConvolvePSFcomputeFFT();
        void StepConvolvePSFcomputeReal();
        void StepConvolvePsfRealRemapArrays(Array<float, 2> &subPixelMap, Array<float, 2> &psfMap);
        void StepConvolvePsfFftRemapArrays(Array<float, 2> &subPixelMap, Array<float, 2> &psfMap);
        void StepConvolvePSFconvolveFFT(Array<float, 2> subPixelMap, Array<float, 2> psfMap, Array<float, 2> &convolvedMap, int threads, bool useFFTWisdom,
		std::string fftWisdomPath);
        void StepConvolvePSFconvolveReal(Array<float, 2> subPixelMap, Array<float, 2> psfMap, Array<float, 2> &convolvedMap);



private:

        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.
        
        std::string convolutionMethod;                       //Parameter retrieved from DataSet.
        int    threads;                                 //Parameter retrieved from DataSet.
        bool   useFFTWisdom;                            //Parameter retrieved from DataSet.
        std::string fftWisdomPath;                           //Parameter retrieved from DataSet.
        int    memoryLimit;                             //Parameter retrieved from DataSet.
        

        Array<float, 2> psfMap;                         //Map array retrieved from DataSet.
        Array<float, 2> pixelMap;                       //Map array retrieved from DataSet.
        Array<float, 2> subPixelMap;                    //Map array retrieved from DataSet.
        Array<float, 2> convolvedMap;                   //Map array for convolving PSF.



};
#endif /* StepConvolvePSF_H_ */