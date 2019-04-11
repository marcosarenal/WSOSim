///////////////////////////////////////////////////////////
// Sky.h: C++ code to define a class that is able to return many properties
//         of the sky and objects in the sky like the sun.
//  Created on:      23-Oct-2012 1:59:58 PM
//  Original author: Joris De Ridder
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


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

#include "Constants.h"
#include "Tabfunction.h"
#include "DataSet.h"
/**
 * Compute zodiacal flux an diffuse stellar background flux for a certain position in the sky using pre-computed tables (skydata.cpp).
 */
class Sky
{
public:

	Sky();
	~Sky();

	double SolarRadiantFlux(const double lambda);
	double SolarRadiantFlux(const double lambda1, const double lambda2);
	double SolarRadiantFlux(vector<double> &lambda, vector<double> &throughput);
	double ZodiacalFlux(const double lambda1, const double lambda2, const double alpha, const double delta);
	double ZodiacalFlux(vector<double> &lambda, vector<double> &throughput, const double alpha, const double delta);
	double StellarBgFlux(const double lambda1, const double lambda2, const double RA, const double decl);
	double StellarBgFlux(vector<double> &lambda, vector<double> &throughput, const double RA, const double decl);

protected:

private:

	const double obliquity; // obliquity of the ecliptic = 23.439 deg. [rad]

	vector<double> integrand;
	Tabfunction<vector<double> > tabfunction;

	void locate(double x, const double *array, int N, int &index);
	void equa2ecl(const double alpha, const double delta, double &lambda, double &beta);

};

