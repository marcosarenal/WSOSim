///////////////////////////////////////////////////////////
//  Constants.cpp
//  Constants parameters
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



#include "Constants.h"

double Constants::e_number       = 2.71828182845904523536029;   // e; Napier's constant or Euler's number 
double Constants::Pi2            = 6.2831853071795862;
double Constants::Pi             = 3.14159265358979323846264;
double Constants::EMPTY          = 99999.12345;
double Constants::pc_c           = 299792458.;                  // speed of light m/s
double Constants::pc_h           = 6.62606876E-34;              // Planck constant
double Constants::pc_k           = 1.3806503E-23;               // Boltzmann constant
double Constants::pc_sigma       = 5.669E-8;                    // Stefan-Boltzmann constant
double Constants::Pi4            = 12.5663706143591725;
double Constants::Pid2           = 1.5707963267948966;
double Constants::Pid4           = 0.7853981633974483;
double Constants::SqrtPi         = 1.7724538509055159;
double Constants::Sqrt2Pi        = 2.5066282746310002;
double Constants::Wuz4Pi         = 3.5449077018110318;
double Constants::Sqrt2          = 1.41421356237309504880169;
double Constants::radiusSun      = 6.9599E5;                    // solar radius in km
double Constants::massSun        = 1.9891E30;                   // solar mass in kg
double Constants::G              = 6.67428E-11;                 // gravitational constant m^3 kg^-1 s^-2
double Constants::RAD2DEG        = 180. / Constants::Pi;        // conversion from radians to degrees
double Constants::DEG2RAD        = Constants::Pi / 180.0;       // conversion from degrees to radians
double Constants::HOUR2RAD       = Constants::Pi / 12.;         // conversion from angular hours to radians
double Constants::DAY2SEC        = 86400.0;                     // conversion from days to seconds
double Constants::PC2M           = 3.0857e16;                   // conversion from parsec to meters
double Constants::M2PC           = 1.0/Constants::PC2M;         // conversion from meters to parsecs
double Constants::RSUN2M         = 6.9599e8;                    // conversion from solar radius to meter
double Constants::MSUN2KG        = 1.989e30;                    // conversion from solar mass to kilogram
double Constants::LSUN2L         = 3.826e26;                    // conversion from solar luminosity to J/s
double Constants::M2RSUN         = 1.0/Constants::RSUN2M;       // conversion from meter to solar radius
double Constants::KG2MSUN        = 1.0/Constants::MSUN2KG;      // conversion from kilograms to solar mass
double Constants::L2LSUN         = 1.0/Constants::LSUN2L;       // conversion from J/s to solar luminosity
double Constants::CLIGHT         = 299792458.0;                 // m/s
double Constants::HPLANCK        = 6.6260755e-34;               // J*s
double Constants::KBOLTZMAN      = 1.380658e-23;                // J/K
double Constants::STEFANBOLTZMAN = 5.67051e-8;                  // J/K^4/m^2/s  (symbol: sigma)
double Constants::AU             = 149597870691.0;              // Astronomical Unit (m)
double Constants::OBLIQUITY      = 0.409087723;                 // obl. of ecliptica (rad)
