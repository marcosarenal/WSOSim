///////////////////////////////////////////////////////////
//  PreProcessingPSF.h
//  Implementation of the Class PreProcessingPSF
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



#ifndef PREPROCESSINGPSF_H_
#define PREPROCESSINGPSF_H_


#include "ParamsPSF.h"
#include "DataSet.h"
#include "LogManager.h"



/**
 * This class performs the duties of the PreProcessing relating to the PSF parameters of any processing mode.
 * this includes calls the ParamsPSF class to set all the PSF parameters in the DataSet.
 */
class PreProcessingPSF
{

public:
        
	PreProcessingPSF();
	PreProcessingPSF(DataSet &m_DataSet);
	virtual ~PreProcessingPSF();
	     
        

private:

        DataSet m_DataSet;
        ParamsPSF m_ParamsPSF;

        
};
#endif /* PREPROCESSINGPSF_H_ */
