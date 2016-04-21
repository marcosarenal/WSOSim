///////////////////////////////////////////////////////////
//  FFT.h
//  FFT methods
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

#ifndef FFT_H
#define FFT_H

#include <fftw3.h>
#include "FileUtilities.h"
#include <sstream>
#include <sys/stat.h>
#include <blitz/array.h>
#include <iostream>

using namespace blitz;
using namespace std;

/**
 * Wrapper class to compute Fast Fourier Transforms using the FFTW++ library.
 */
class FFT
{

public:
	/**
	 * Empty constructor. Does nothing.
	 */
	FFT();

	/**
	 * 1D FFT transformation forward.
	 */
	static void forward(Array<float, 1> &in, Array<complex<float> , 1> &out);

	/**
	 * 1D FFT transformation backward.
	 */
	static void backward(Array<complex<float> , 1> &in, Array<float, 1> &out);//1D FFT transformation backward

	/**
	 * 2D FFT transformation forward.
	 * \param in 2D input array
	 * \param out Transformed 2D output array
	 * \param threads Number of treads that should be used in parallel for the computation of the FFT. Does not work well.
	 * \param useFFTWisdom Set to true when FFTW-wisdom should be used. Wisdom consists of a library of pre-computed matrices that speed up computation of the FFT.
	 * \param fftWisdomPath Path on disk where the wisdom information is stored.
	 */
	static void forward(Array<float, 2> &in, Array<complex<float> , 2> &out, int threads, bool useFFTWisdom = false, string fftWisdomPath = "");

	/**
	 * 2D FFT transformation backward.
	 * \param colst2 Number of columns of real 2D array. Differs from ncol of complex input array (see FFTW manual) http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
	 * \param in 2D input array
	 * \param out Transformed 2D output array
	 * \param threads Number of treads that should be used in parallel for the computation of the FFT. Does not work well.
	 * \param useFFTWisdom Set to true when FFTW-wisdom should be used. Wisdom consists of a library of pre-computed matrices that speed up computation of the FFT.
	 * \param fftWisdomPath Path on disk where the wisdom information is stored.
	 */
	static void backward(int colst2, Array<complex<float> , 2> &in, Array<float, 2> &out, int threads, bool useFFTWisdom = false,
			string fftWisdomPath = "");

private:

};
#endif
