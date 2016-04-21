///////////////////////////////////////////////////////////
//  Controller.h
//  Implementation of the Class Controller
//  Created on:      23-Oct-2012 1:59:57 PM
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



#ifndef CONTROLLER_H_
#define CONTROLLER_H_

#include <string>
#include <iostream>

#include "GlobalVariables.h"
#include "LogManager.h"
#include "PreProcessingCommon.h"
#include "PreProcessingCCD.h"
#include "PreProcessingWUVS.h"
#include "PreProcessingPSF.h"
#include "PreProcessingCMOS.h"
#include "PreProcessingPhotometry.h"
#include "ProcessingCCD.h"
#include "ProcessingWUVS.h"
#include "ProcessingCMOS.h"
#include "ProcessingPhotometry.h"
#include "PostProcessing.h"
#include "DataSet.h"

using namespace std;


/**
 * This class is in charge of controlling and monitorize the whole system. Here is
 * selected the type of processing to be performed depending on what is the
 * simulation selected by the user (CCD, CMOS, photometry). According to this
 * selection, different classes are triggered in the PreProcessing, Processing and
 * PostProcessing subsystems. When each of the calls are performed, the status is
 * reported to the LogManager.
 */
class Controller
{

public:	
    
    Controller();
    virtual ~Controller();
        
    static void runController(int argc, char ** argv);
    static void runCCDandPhotometryController(string parameterFile, string photometryParameterFile);
    static void runCMOSandPhotometryController(string parameterFile, string photometryParameterFile);
    static void runCCDController(string parameterFile);
    static void runCMOSController(string parameterFile);
    static void runPhotometryController(string parameterFile, string photometryParameterFile);
    static void runWUVSController(string parameterFile);

private:
        


};
#endif /* CONTROLLER_H_ */