///////////////////////////////////////////////////////////
//  LogManager.h
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




#ifndef LOGMANAGER_H_
#define LOGMANAGER_H_

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

/**
 * Performs and controls the logs of the whole system.
 */
class LogManager
{
public:
	LogManager()
	{
		showLogsFlag = true;
	};
        
	LogManager(const LogManager& temp);
	LogManager &operator=(const LogManager &temp);
	void LogManagerGenerateLogFile(std::string fileName);

	friend stringstream &operator<<(stringstream & out, const std::string & outp);

	void LogManagerSetFile(std::string fileName);
	void remove();
        
        
	void LogManagerAppendLog(stringstream &text);
	void LogManagerAppendLog()
	{
		LogManagerAppendLog(log);
		log.str("");
	};
        
        
	void LogManagerAppendLogAndShow(stringstream &text);
	void LogManagerAppendLogAndShow()
	{
		LogManagerAppendLogAndShow(log);
		log.str("");
	};
        
        
	void LogManagerShowLog(stringstream &text);
	void LogManagerShowLog()
	{
		LogManagerShowLog(log);
		log.str("");
	};

	void showLogs(bool b)
	{
		showLogsFlag = b;
	};

	static void LogManagerShowLog(std::string text);
	static bool showLogsFlag;
	static stringstream log;
    void LogManagerPresentation();

      
private:
	std::string fileName;
};
#endif /* LOGMANAGER_H_ */