///////////////////////////////////////////////////////////
//  Statistics.cpp
//  Statistical methods
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



#include "Statistics.h"
#include <cmath>
#include "MathTools.h"

Statistics::Statistics()
{
	// TODO
}

//==============================================================================




//==============================================================================
/**
 * @param x
 * @param mean
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x, double mean)
{
	double dev = 0.;

	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += pow(mean - x(i), 2);
	}
	dev = sqrt(dev / (double(x.size()) - 1.));

	return dev;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param w
 * @param mean
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x, Array<double, 1> w, double mean)
{
	double dev = 0., wsum = 0.;

	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += (pow(mean - x(i), 2) * w(i));
		wsum += w(i);
	}
	dev = sqrt(dev / wsum);

	return dev;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param w
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x, Array<double, 1> w)
{
	double dev = 0., wsum = 0.;
	double meanv = mean(x);
	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += pow2(meanv - x(i));
		wsum += w(i);
	}

	dev = sqrt(dev / wsum);
	return dev;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param mean
 * @param indices
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x, double mean, Array<int, 1> indices)
{
	double dev = 0.;

	for (unsigned int i = 0; i < indices.size(); i++)
	{
		dev += pow(mean - x(indices(i)), 2);
	}
	dev = sqrt(dev / (double(indices.size()) - 1.));

	return dev;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @return 
 */
double Statistics::getMean(Array<double, 1> x)
{
	double ave = 0.;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		ave += x(i);
	}
	ave /= x.size();
	return ave;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param w
 * @return 
 */
double Statistics::getMean(Array<double, 1> x, Array<double, 1> w)
{
	double ave = 0., wSum = 0.;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		ave += x(i) * w(i);
		wSum += w(i);
	}
	ave /= wSum;
	return ave;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param indices
 * @return 
 */
double Statistics::getMean(Array<double, 1> x, Array<int, 1> indices)
{
	double ave = 0.;
	for (unsigned int i = 0; i < indices.size(); i++)
	{
		ave += x(indices(i));
	}
	ave /= indices.size();
	return ave;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param xMean
 * @param yMean
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x, Array<double, 1> y, double xMean, double yMean)
{
	double dev = 0.;

	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += (pow(x(i) - xMean, 2) + pow(y(i) - yMean, 2));
	}
	dev = sqrt(dev / (double(x.size()) - 1.));

	return dev;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param x0
 * @param y0
 * @return 
 */
double Statistics::distance(double x, double y, double x0, double y0)
{
	return (sqrt(pow(x - x0, 2.) + pow(y - y0, 2.)));
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @return 
 */
double Statistics::getMin(Array<double, 1> x)
{
	double min = x(0);
	for (unsigned int i = 1; i < x.size(); i++)
		if (x(i) < min)
			min = x(i);

	return min;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @return 
 */
double Statistics::getMax(Array<double, 1> x)
{

	double max = x(0);
	for (unsigned int i = 1; i < x.size(); i++)
		if (x(i) > max)
			max = x(i);

	return max;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param mean
 * @param seedRNG
 * @return 
 */
double Statistics::getPoisson(double mean, unsigned int seedRNG)
{
	if (mean <= 0)
		return 0;

	else if (mean < 30) //poisson
	{
		UniformClosed<double, ranlib::MersenneTwister, ranlib::independentState> rand;
		rand.seed(seedRNG);

		double L = exp(-mean);
		double p = 1;
		int k = 0;

		do
		{
			k += 1;
			p *= rand.random();
		} while (p >= L);

		return (k - 1);
	}

	else //normal
	{
		Normal<double> rand(mean, sqrt(mean));
		rand.seed(seedRNG);

		int number = int(round(rand.random()));
		if (number < 0)
			number = 0;
		return (number);

	}
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @return 
 */
double Statistics::getStdDev(Array<double, 1> x)
{
	double dev = 0.;
	double meanv = mean(x);
	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += pow2(meanv - x(i));
	}
	dev = sqrt(dev / (double(x.size()) - 1.));

	return dev;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @return 
 */
double Statistics::getStdDev(Array<float, 1> x)
{
	float dev = 0.;
	float meanv = mean(x);
	for (unsigned int i = 0; i < x.size(); i++)
	{
		dev += pow2(meanv - x(i));
	}
	dev = sqrt(dev / (float(x.size()) - 1.));

	return dev;

}
//==============================================================================

