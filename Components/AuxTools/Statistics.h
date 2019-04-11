///////////////////////////////////////////////////////////
//  Statistics.h
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



#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "random/uniform.h"
#include "random/normal.h"
#include "blitz/array.h"

using namespace ranlib;
using namespace blitz;

class Statistics
{

public:
	Statistics();
	static double getStdDev(Array<double, 1> x);
	static double getStdDev(Array<float, 1> x);
	static double getStdDev(Array<double, 1> x, double mean, Array<int, 1> indices);
	static double getStdDev(Array<double, 1> x, double mean);
	static double getStdDev(Array<double, 1> x, Array<double, 1> w);
	static double getStdDev(Array<double, 1> x, Array<double, 1> w, double mean);
	static double getMean(Array<double, 1> x, Array<int, 1> indices);
	static double getMean(Array<double, 1> x);
	static double getMean(Array<double, 1> x, Array<double, 1> w);
	static double getStdDev(Array<double, 1> x, Array<double, 1> y, double xMean, double yMean);
	static double distance(double x, double y, double x0, double y0);
	static double getMin(Array<double, 1> x);
	static double getMax(Array<double, 1> x);
	static double getPoisson(double mean, unsigned int seedRNG); //return a poisson distributed value around mean. for values larger that 30 we take a normal distribution
	//  static double getStdDev ( Array<double, 1> x, Array<double, 1> w );
};
#endif /* STATISTICS_H_ */
