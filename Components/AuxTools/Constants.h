///////////////////////////////////////////////////////////
//  Constants.h
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



#ifndef CONSTANTS_H_
#define CONSTANTS_H_


class Constants  
{

public:
	Constants() {};

        static double Pi2;
        static double Pi;
        static double pc_c;                     // speed of light m/s
        static double pc_h;                     // Planck constant
        static double pc_k;                     // Boltzmann constant
        static double pc_sigma;                 // Stefan-Boltzmann constant
        static double Pi4;
        static double Pid2;
        static double Pid4;
        static double SqrtPi;
        static double Sqrt2Pi;
        static double Wuz4Pi;
        static double Sqrt2;
        static double radiusSun;                // solar radius in km
        static double massSun;                  // solar mass in kg
        static double G;                        // gravitational constant m^3 kg^-1 s^-2
	static double RAD2DEG;                  // conversion from radians to degrees
        static double DEG2RAD;                  // conversion from degrees to radians
        static double HOUR2RAD;                 // conversion from angular hours to radians
        static double DAY2SEC;                  // conversion from days to seconds
        static double PC2M;                     // conversion from parsec to meters
        static double M2PC;                     // conversion from meters to parsecs
        static double RSUN2M;                   // conversion from solar radius to meter
        static double MSUN2KG;                  // conversion from solar mass to kilogram
        static double LSUN2L;                   // conversion from solar luminosity to J/s
        static double M2RSUN;                   // conversion from meter to solar radius
        static double KG2MSUN;                  // conversion from kilograms to solar mass
        static double L2LSUN;                   // conversion from J/s to solar luminosity
        static double CLIGHT;                   // m/s
        static double HPLANCK;                  // J*s
        static double KBOLTZMAN;                // J/K
        static double STEFANBOLTZMAN;           // J/K^4/m^2/s  (symbol: sigma)
        static double AU;                       // Astronomical Unit (m)
        static double OBLIQUITY;                // obl. of ecliptica (rad)
	
        static double EMPTY;

};
#endif /* CONSTANTS_H_ */
