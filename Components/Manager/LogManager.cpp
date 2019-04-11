///////////////////////////////////////////////////////////
//  LogManager.cpp
//  Implementation of the Class LogManager
//  Created on:      25-Nov-2010 
//  Original author: zima
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



#include "LogManager.h"


std::stringstream LogManager::log;
bool LogManager::showLogsFlag = true;


//==============================================================================
/**
 * Constructor method
 */

//==============================================================================
/**
 * 
 * @param temp
 */
LogManager::LogManager(const LogManager& temp)
{
	fileName = temp.fileName;
	showLogsFlag = temp.showLogsFlag;
}
//==============================================================================

//==============================================================================
/**
 *  Copy assignment constructor
 */
LogManager &LogManager::operator=(const LogManager &temp)
{

	if (this != &temp)
	{
		fileName = temp.fileName;
		showLogsFlag = temp.showLogsFlag;
	}

	return *this;
}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================

/**
 * This function generates a log file to keep a log of all the processing information
 * generated during the simulation.
 * 
 * @param fileName File name including the absolute path of the log file to be created.
 */
void LogManager::LogManagerGenerateLogFile(std::string fileName)
{
	showLogsFlag = true;
	LogManagerSetFile(fileName);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 */
void LogManager::LogManagerSetFile(std::string fileName)
{
//        LogManager::log << "        test 1 " ;

	this->fileName = fileName;
//                LogManager::log << "        test 2 " ;

	remove();
//                LogManager::log << "        test 3 " ;

}
//==============================================================================




//==============================================================================
/**
 * 
 */
void LogManager::remove()
{
//	unlink(fileName.c_str());
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param text
 */
void LogManager::LogManagerAppendLog(std::stringstream &text)
{
//	if (GlobalVariables::rankMPI > 0)
//		return;

	if (!LogManager::showLogsFlag)
		return;

	ofstream out(fileName.c_str(), ios_base::out | ios_base::app);
	if (!out.is_open())
	{
		std::cerr << "\nError (LogManager::LogManagerAppendLog()): Unable to open output log file " << fileName;
		exit(1);
	}
	out.precision(11);
	out << text.str() << std::endl;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param text
 */
void LogManager::LogManagerShowLog(std::stringstream &text)
{
//	if (GlobalVariables::rankMPI > 0)
//		return;

	if (LogManager::showLogsFlag)
		std::cout << text.str() << std::endl;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param log
 */
void LogManager::LogManagerShowLog(std::string text)
{
//	if (GlobalVariables::rankMPI > 0)
//		return;

	if (LogManager::showLogsFlag)
		std::cout << text << std::endl;
}
//==============================================================================





//==============================================================================
/**
 * 
 * @param text
 */
void LogManager::LogManagerAppendLogAndShow(std::stringstream &text)
{
	LogManagerAppendLog(text);
	LogManagerShowLog(text);
}
//==============================================================================



//==============================================================================
/**
 * Credits for the presentation of the program. 
 */
void LogManager::LogManagerPresentation()
{
    log << "     "<<std::endl;
    log << "     *****************************************************************************"<<std::endl;
    log << "     *                                  WSOSim                                   *"<<std::endl;
    log << "     *                          WSO End-to-End Simulator                         *"<<std::endl;
    log << "     *                               Version 3.0                                 *"<<std::endl;
    log << "     *                                11/05/2016                                 *"<<std::endl;
    log << "     *                         Based on the PLATO Simulator                      *"<<std::endl;
    log << "     *                                                                           *"<<std::endl;
    log << "     *                        Developer: Pablo Marcos Arenal                     *"<<std::endl;
    log << "     *                                    UCM                                    *"<<std::endl;
    log << "     *                               pablmar@ucm.es                              *"<<std::endl;
    log << "     *****************************************************************************"<<std::endl;	
    log << "     "<<std::endl;

}
//==============================================================================

                                                                       