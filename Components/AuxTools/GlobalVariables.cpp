///////////////////////////////////////////////////////////
//  GlobalVariables.cpp
//  Implementation of Global Variables
//  Created on:      14-Mar-2013 1:59:58 PM
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




#include "GlobalVariables.h"








//==============================================================================
/**
 * Constructor method
 */
GlobalVariables::GlobalVariables(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
GlobalVariables::~GlobalVariables(){}
//==============================================================================




//Initialization of the LogManager.        
LogManager GlobalVariables::logManager = LogManager();

