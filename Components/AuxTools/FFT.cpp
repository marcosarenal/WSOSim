///////////////////////////////////////////////////////////
//  FFT.cpp
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


#include "FFT.h"


FFT::FFT()
{
	// TODO
}
//


//****************************************************************************80
void FFT::forward(Array<float, 1> &in, Array<complex<float> , 1> &out)
{
	int n = in.size();
	fftw_complex *fin, *fout;
	fftw_plan p;

	fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

	for (int i = 0; i < n; i++)
	{
		fin[i][0] = in(i);
		fin[i][1] = 0.;
	}

	p = fftw_plan_dft_1d(n, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE );

	fftw_execute(p);

	out.resize(n);
    out = 0.0;

	for (int i = 0; i < n; i++)
	{
		out(i).real() = fout[i][0];
		out(i).imag() = fout[i][1];
	}

	fftw_destroy_plan(p);
	fftw_free(fin);
	fftw_free(fout);

}

//****************************************************************************80
void FFT::backward(Array<complex<float> , 1> &in, Array<float, 1> &out)
{
	int n = in.size();
	fftw_complex *fin, *fout;
	fftw_plan p;

	fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

	for (int i = 0; i < n; i++)
	{
		fin[i][0] = in(i).real();
		fin[i][1] = in(i).imag();
	}

	p = fftw_plan_dft_1d(n, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE );

	fftw_execute(p);

	out.resize(n);
    out = 0.0;

	for (int i = 0; i < n; i++)
		out(i) = fout[i][0];

	fftw_destroy_plan(p);
	fftw_free(fin);
	fftw_free(fout);

}

//****************************************************************************
void FFT::forward(Array<float, 2> &in, Array<complex<float> , 2> &out, int threads, bool useFFTWisdom, std::string fftWisdomPath)
{
	//qDebug ( "void FFT::forward ( Array<float, 2> in, Array<complex<float>, 2 > &out )" );
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	//  fftw_set_timelimit ( 300 );
	int rows = in.rows();
	int cols = in.cols();
	int colsd2 = cols / 2 + 1;

	if (!FileUtilities::dirExists(fftWisdomPath.c_str()))
	{
		std::cerr << "\nError (FFT::forward): Wisdom directory '" << fftWisdomPath << "' does not exist. Please create it." << std::endl;
        exit(1);
	}

	std::stringstream ss;
	ss << fftWisdomPath << "/fftwisdomr2c_r" << rows << "_c" << cols << ".dat";

	std::string wisdomFile = ss.str();

	double * fin = (double*) fftw_malloc(sizeof(double) * rows * cols);
	fftw_complex * fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * colsd2);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			fin[cols * i + j] = in(i, j);

	in.free();
	fftw_plan p;

	if (useFFTWisdom)
	{
		if (FileUtilities::fileExists(wisdomFile))
		{
			FILE * pFile;
			pFile = fopen(wisdomFile.c_str(), "r");
			fftw_import_wisdom_from_file(pFile);
			fclose(pFile);

			p = fftw_plan_dft_r2c_2d(rows, cols, fin, fout, FFTW_WISDOM_ONLY );

			if (p == NULL)
			{
	
                LogManager::log<< "Creating a new r2c wisdom plan for FFT convolution (might take some time).";
                LogManager::log<< "Once the plan is available, subsequent FFTs using this plan perform very fast." << std::endl;
                GlobalVariables::logManager.LogManagerShowLog();  

				p = fftw_plan_dft_r2c_2d(rows, cols, fin, fout, FFTW_MEASURE );
				pFile = fopen(wisdomFile.c_str(), "w");
				fftw_export_wisdom_to_file(pFile);
				fclose(pFile);
			}
		}
		else
		{
                
            LogManager::log<< "Creating new FFTW wisdom r2c file."<<std::endl;
            LogManager::log<< "Computing new r2c wisdom plan for FFT convolution (might take some time - for large image dimensions perhaps get a coffee). ";
            LogManager::log<< "Once the plan is available, subsequent FFTs using this plan perform very fast. "<<std::endl;
            LogManager::log<< ""<<std::endl;
            LogManager::log<< "Another option would be to turn off the usage of wisdom by setting <UseFFTWisdom>0</UseFFTWisdom>  ";
            LogManager::log<< "in the parameter file, or use real space convolution (<Convolution>Real</Convolution>).";
            GlobalVariables::logManager.LogManagerShowLog(); 
            
            
            
			p = fftw_plan_dft_r2c_2d(rows, cols, fin, fout, FFTW_MEASURE );

			FILE * pFile;
			pFile = fopen(wisdomFile.c_str(), "w");
			fftw_export_wisdom_to_file(pFile);
			fclose(pFile);
		}
	}
	else
		p = fftw_plan_dft_r2c_2d(rows, cols, fin, fout, FFTW_ESTIMATE );

	fftw_execute(p);

	fftw_free(fin);

	out.resize(rows, colsd2);
    out = 0.0;


	for (int i = 0; i < rows; i++)
		for (int j = 0; j < colsd2; j++)
		{
			out(i, j).real() = fout[colsd2 * i + j][0];
			out(i, j).imag() = fout[colsd2 * i + j][1];
		}

	fftw_destroy_plan(p);
	fftw_free(fout);
	fftw_cleanup_threads();
}

//****************************************************************************80
void FFT::backward(int colst2, Array<complex<float> , 2> &in, Array<float, 2> &out, int threads, bool useFFTWisdom, std::string fftWisdomPath)
//colst2: number of columns of real 2D array. Differs from ncol of complex input array (see FFTW manual)
//http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
{
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	//  fftw_set_timelimit ( 300 );
	int rows = in.rows();
	int cols = in.cols();

	if (!FileUtilities::dirExists(fftWisdomPath.c_str()))
	{
		std::cerr << "\nError (FFT::backward): Wisdom directory " << fftWisdomPath << " does not exist." << std::endl;
		exit(1);
	}

	std::stringstream ss;
	ss << fftWisdomPath << "/fftwisdomc2r_r" << rows << "_c" << colst2 << ".dat";
	std::string wisdomFile = ss.str();

	fftw_complex * fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * cols);
	double * fout = (double*) fftw_malloc(sizeof(double) * rows * colst2);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			fin[cols * i + j][0] = in(i, j).real();
			fin[cols * i + j][1] = in(i, j).imag();
		}

	in.free();
	fftw_plan p;
	if (useFFTWisdom)
	{
		if (FileUtilities::fileExists(wisdomFile))
		{
			FILE * pFile;
			pFile = fopen(wisdomFile.c_str(), "r");
			fftw_import_wisdom_from_file(pFile);
			fclose(pFile);

			p = fftw_plan_dft_c2r_2d(rows, colst2, fin, fout, FFTW_WISDOM_ONLY );

			if (p == NULL)
			{
                LogManager::log << "Computing new c2r wisdom plan for FFT convolution (might take some time).";
                LogManager::log << "Once the plan is available, subsequent FFTs using this plan perform very fast.";        
                GlobalVariables::logManager.LogManagerShowLog(); 

				p = fftw_plan_dft_c2r_2d(rows, colst2, fin, fout, FFTW_MEASURE );
				pFile = fopen(wisdomFile.c_str(), "w");
				fftw_export_wisdom_to_file(pFile);
				fclose(pFile);
			}
		}
		else
		{
            LogManager::log << "    Creating new FFTW wisdom c2r file." <<std::endl;
            LogManager::log << "    Computing new c2r wisdom plan for FFT convolution (might take some time - for large image dimensions perhaps get a numPrescanRows)."<<std::endl;
            LogManager::log << "    Once the plan is available, subsequent FFTs using this plan perform very fast."<<std::endl;
            LogManager::log << "    Another option would be to turn off the usage of wisdom by setting <UseFFTWisdom>0</UseFFTWisdom> in the parameter file, "<<std::endl;
            LogManager::log << "    or use real space convolution (<Convolution>Real</Convolution>).";
            GlobalVariables::logManager.LogManagerShowLog(); 

			p = fftw_plan_dft_c2r_2d(rows, colst2, fin, fout, FFTW_MEASURE );

			FILE * pFile;
			pFile = fopen(wisdomFile.c_str(), "w");
			fftw_export_wisdom_to_file(pFile);
			fclose(pFile);
		}
	}
	else
		p = fftw_plan_dft_c2r_2d(rows, colst2, fin, fout, FFTW_ESTIMATE );

	fftw_execute(p);

	out.resize(rows, colst2);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < colst2; j++)
			out(i, j) = fout[colst2 * i + j];

	fftw_destroy_plan(p);
	fftw_free(fin);
	fftw_free(fout);
	fftw_cleanup_threads();

}

