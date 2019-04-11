///////////////////////////////////////////////////////////
//  ParamsPhotometry.h
//  Implementation of the Class ParamsPhotometry
//  Created on:      23-Oct-2012 1:59:59 PM
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



#ifndef PARAMSPHOTOMETRY_H_
#define PARAMSPHOTOMETRY_H_

#include "DataSet.h"

class ParamsPhotometry
{

public:
	DataSet* m_DataSet;

	ParamsPhotometry();
	virtual ~ParamsPhotometry();
	void paramsPhotometryCalculation();

};
#endif /* PARAMSPHOTOMETRY_H_ */
