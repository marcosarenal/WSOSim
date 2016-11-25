///////////////////////////////////////////////////////////
//  FileUtilities.h
//  File utilities methods
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


#ifndef FILEUTILITIES_H_
#define FILEUTILITIES_H_

#include <iostream>
#include <string>
#include <vector>
#include "DataSet.h"
#include "CCfits/CCfits"
#include "blitz/array.h"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include <glob.h>
#include <dirent.h>

using namespace std;
using namespace CCfits;


/**
 * Convenience functions to read and write files consisting of several columns.
 */
class FileUtilities
{

public:
	FileUtilities()
	{
	}
	;
	static int countFiles(string pattern);
	static string readAll(string fileName);
	static bool readCol(string fileName, int col, vector<string> &var, int skip = 0);
	static bool readCol(string fileName, int col, vector<double> &var, int skip = 0);
	static bool readCol(string fileName, int col, vector<float> &var, int skip = 0);
	static bool readCol(string fileName, int col, Array<double, 1> &var, int skip = 0);
	static bool readCol(string fileName, int col, Array<float, 1> &var, int skip = 0);
	static double readValue(string file, int column, int row);
	static int countLines(string fileName);

	template<class T>
	static bool from_string(T& t, const string& s, std::ios_base& (*f)(std::ios_base&))
	{
		std::istringstream iss(s);
		return !(iss >> f >> t).fail();
	}
	static bool doMkdir(const char *path, mode_t mode);
	static bool dirExists(const char *path);
	static bool fileExists(string strFilename);
	static void getDir(string dir, vector<string> &files, string filter = "");
	static string getCompleteBaseName(const string& str);
	static string getCompletePath(const string& str);
	static string removeExtension(const string& str);
        
    static void FileUtilitiesCutAndWriteFITS(double exposureTime, int cut, int numPrescanRows, int numSmearingOverscanRows, Array<float, 2>  biasRegisterMap, 
                                             Array<float, 2>  smearingMap, Array<float, 2> &map, string fileName);
    static void FileUtilitiesWriteFITS(double exposureTime, int cut, Array<double, 2> &map, string fileName);
    static void FileUtilitiesWriteFITS(double exposureTime, int cut, Array<float, 2> &map, string fileName);
   
    static void FileUtilitiesWriteFits( string fileName, Array<double, 2> map, int cut = 0 );
    static void FileUtilitiesWriteFits( string fileName, Array<float, 2> map, int cut = 0 );   
    
    static void FileUtilitiesReadFits( string fileName, Array<int, 2> &flux );
    static void FileUtilitiesReadFits( string fileName, Array<double, 2> &flux );
    static void FileUtilitiesReadFits( string fileName, Array<float, 2> &flux );

    static void FileUtilitiesReadExternalFits( string fileName, Array<float, 2> &flux );


};
#endif /* FILEUTILITIES_H_ */
