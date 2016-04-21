///////////////////////////////////////////////////////////
//  ParamsCMOS.h
//  Implementation of the Class ParamsCMOS
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



#ifndef PARAMSCMOS_H_
#define PARAMSCMOS_H_

#include "DataSet.h"

/**
 * Defines all properties of the CMOS, contains the sub- and normal-pixel maps and
 * contains several methods and algorithms for dealing with a CCD. Describes the
 * physical and electrical properties of the CMOS. In the simulations, it is
 * assumed that the CMOS of all telescopes have the same general properties. The
 * read-out direction of the CMOSis assumed to be oriented in negative y-direction
 * and the read-out strip is below the y=0 row. The field of view of the
 * CMOSdetermines which stars affect the sub-field through read-out smearing (for
 * a charge transfer time > 0 s). The flatfield ( pixel non-uniform response) is
 * computed by considering a spatial 1/f-response of the sensitivity. Due to the
 * consideration of computing time, not the complete CMOS is modelled with sub-
 * pixel precision but only a small sub-field with a side length of a few hundred
 * pixels.
 */
class ParamsCMOS
{

public:
	DataSet* m_DataSet;

	ParamsCMOS();
	virtual ~ParamsCMOS();
	void paramsCMOSCalculation();

};
#endif /* PARAMSCMOS_H_ */
