///////////////////////////////////////////////////////////
//  FileUtilities.cpp
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




#include "FileUtilities.h"




//==============================================================================
//Functions:
//==============================================================================

/**
 * 
 * @param pattern
 * @return 
 */
int FileUtilities::countFiles ( std::string pattern )
{
   glob_t globbuf;
   globbuf.gl_offs = 2;

   return glob ( pattern.c_str(), GLOB_DOOFFS, NULL, &globbuf );

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @return 
 */
int FileUtilities::countLines ( std::string fileName )
{
   int counter = 0;
   std::string line;

   ifstream myfile ( fileName.c_str() );
   if ( myfile.is_open() )
   {
      while ( getline ( myfile, line ) )
      {
         counter++;
      }

      myfile.close();
      return counter;
   }
   else
   {
       LogManager::log << "Unable to open file " << fileName << std::endl;
       GlobalVariables::logManager.LogManagerShowLog(); 
       return 0;
   }
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @return 
 */
std::string FileUtilities::readAll ( std::string fileName )
{
   std::string fileContent, line;
   ifstream myfile ( fileName.c_str() );
   if ( myfile.is_open() )
   {
      while ( getline ( myfile, line ) )
      {
         fileContent += line + "\n";
      }
      myfile.close();
      return fileContent;
   }
   else
   {
       std::cerr << "\nError (FileUtilities::readAll()): Unable to open file " << fileName << std::endl;
       exit (1);
   }
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @param col
 * @param var
 * @param skip
 * @return 
 */
bool FileUtilities::readCol ( std::string fileName, int col, vector<std::string> &var, int skip )
{
   var.clear();
   string str, line;
   istringstream iss;
   ifstream myfile ( fileName.c_str() );
   if ( myfile.is_open() )
   {
      while ( getline ( myfile, line ) )
      {
         iss.clear();
         iss.str ( line );
         for ( int i = 0;i <= col;i++ )
            iss >> str;
         var.push_back ( str );
      }
      myfile.close();
      return true;
   }
   else
   {
       std::cerr << "\nError (FileUtilities::readCol()): Unable to open file " << fileName << std::endl;
       exit (1);
   }

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @param col
 * @param var
 * @param skip
 * @return 
 */
bool FileUtilities::readCol ( std::string fileName, int col, vector<double> &var, int skip )
{
   var.clear();
   std::string str, line;
   double value;

   istringstream iss;
   ifstream myfile ( fileName.c_str() );
   if ( myfile.is_open() )
   {
      while ( getline ( myfile, line ) )
      {
         iss.clear();
         iss.str ( line );
         for ( int i = 0;i <= col;i++ )
            iss >> str;

         from_string<double> ( value, str, std::dec );
         var.push_back ( value );
      }
      myfile.close();
      return true;
   }
   else
   {
       std::cerr << "\nError (FileUtilities::readCol()): Unable to open file " << fileName << std::endl;
       exit (1);
   }


}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @param col
 * @param var
 * @param skip
 * @return 
 */
bool FileUtilities::readCol ( std::string fileName, int col, vector<float> &var, int skip )
{
   var.clear();
   std::string str, line;
   double value;

   istringstream iss;
   ifstream myfile ( fileName.c_str() );
   if ( myfile.is_open() )
   {
      while ( getline ( myfile, line ) )
      {
         iss.clear();
         iss.str ( line );
         for ( int i = 0;i <= col;i++ )
            iss >> str;

         from_string<double> ( value, str, std::dec );
         var.push_back ( value );
      }
      myfile.close();
      return true;
   }
   else
   {
       std::cerr << "\nError (FileUtilities::readCol()): Unable to open file " << fileName << std::endl;
       exit (1);
   }


}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @param col
 * @param var
 * @param skip
 * @return 
 */
bool FileUtilities::readCol ( std::string fileName, int col, Array<double, 1> &var, int skip )
{
   vector<double> vec;
   readCol ( fileName, col, vec, skip );
   var.resize ( vec.size() );
   var = 0.0;
   
   for ( unsigned int i = 0; i < vec.size(); i++ )
   {
      var ( i ) = vec[i];
   }
    return true;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fileName
 * @param col
 * @param var
 * @param skip
 * @return 
 */
bool FileUtilities::readCol ( std::string fileName, int col, Array<float, 1> &var, int skip )
{
   vector<float> vec;
   readCol ( fileName, col, vec, skip );
   var.resize ( vec.size() );
   var = 0.0;

   for ( unsigned int i = 0; i < vec.size(); i++ )
   {
       var ( i ) = vec[i];
   }
   return true;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param path
 * @param mode
 * @return 
 */
bool FileUtilities::doMkdir ( const char *path, mode_t mode )
{
   struct stat st;
   bool status = true;

   if ( stat ( path, &st ) != 0 )
   {
      /* Directory does not exist */
      if ( mkdir ( path, mode ) != 0 )
         status = false;
   }
   else if ( !S_ISDIR ( st.st_mode ) )
   {
      errno = ENOTDIR;
      status = false;
   }

   return ( status );
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param path
 * @return 
 */
bool FileUtilities::dirExists ( const char *path )
{
   struct stat st;

   if ( stat ( path, &st ) != 0 )
      return false;

   return true;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param strFilename
 * @return 
 */
bool FileUtilities::fileExists ( std::string strFilename )
{
   struct stat stFileInfo;

   // Attempt to get the file attributes
   if ( stat ( strFilename.c_str(), &stFileInfo ) == 0 )
   {
      // We were able to get the file attributes
      // so the file obviously exists.
      return true;
   } else
   {
      // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.
      return false;
   }

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param file
 * @param column
 * @param row
 * @return 
 */
double FileUtilities::readValue ( std::string file, int column, int row )
{

   std::string str, line;
   double value;
   istringstream iss;
   ifstream myfile ( file.c_str() );
   if ( myfile.is_open() )
   {
      int r = 0;
      while ( getline ( myfile, line ) && r < row )
         r++;
      iss.clear();
      iss.str ( line );
      for ( int i = 0;i <= column;i++ )
         iss >> str;
      from_string<double> ( value, str, std::dec );
      myfile.close();
      return value;
   }
   else
   {
      std::cerr << "\nError (FileUtilities::readValue()): Unable to open file " << file << std::endl;
      exit ( 1 );
   }

}
//==============================================================================




//==============================================================================
/**
 * This function returns the list of files in a directory
 * @param dir 
 * @param files
 * @param filter
 */
void FileUtilities::getDir ( std::string dir, vector<std::string> &files, std::string filter )
{
   DIR *dp;
   std::string fileName;
   struct dirent *dirp;
   if ( ( dp  = opendir ( dir.c_str() ) ) == NULL ) {
      std::cerr << "\nError (" << errno << ") (FileUtilities::getDir) opening " << dir << std::endl;
      exit ( 1 );
   }

   while ( ( dirp = readdir ( dp ) ) != NULL )
   {
      fileName = std::string ( dirp->d_name );
      if ( filter == "" )
         files.push_back ( fileName );
      else if ( fileName.find ( filter ) != std::string::npos )
         files.push_back ( fileName );

   }
   closedir ( dp );
}
//==============================================================================




//==============================================================================
/**
 * Get the file name (e.g. if complete path is "/home/user/simulations/test" then this will return the std::string "test"
 * @param str
 * @return 
 */
std::string FileUtilities::getCompleteBaseName ( const std::string& str )
{
   size_t found;
   found = str.find_last_of ( "/\\" );
   return str.substr ( found + 1 );
}
//==============================================================================




//==============================================================================
/**
 * Get the directory name (e.g. if complete path is "/home/user/simulations/test" then this will return the std::string "/home/user/simulations/"
 * @param str
 * @return 
 */
std::string FileUtilities::getCompletePath ( const std::string& str )
{
   size_t found;
   found = str.find_last_of ( "/\\" );
   return str.substr ( 0, found );
}
//==============================================================================




//==============================================================================
/**
 * Get the file name (e.g. if complete path is "/home/user/simulations/test.txt" then this will return the std::string "/home/user/simulations/test"
 * @param str
 * @return 
 */
std::string FileUtilities::removeExtension ( const std::string& str )
{
   size_t found;
   found = str.find_last_of ( "." );
   return str.substr ( 0, found );
}
//==============================================================================




//==============================================================================
/**
 * Wrapper class to handle fits files. Provides convenience methods to write fits files from a Blitz array of floats.
 * 
 * @param exposureTime Exposure time, to be written in the FITS file header.
 * @param cut Controls how much of the edge of the image is not written to FITS file (should be half size of PSF).
 * @param &map Blitz Array<float, 2> to be written in the FITS file.
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 */
void FileUtilities::FileUtilitiesWriteFITS(double exposureTime, int cut, Array<float, 2> &map, std::string fileName)
{
         
        LogManager::log << "    Creating FITS file " ;
        GlobalVariables::logManager.LogManagerShowLog();        

        
        fileName += ".fits";
        long naxis = 2;
        int cutx = cut;
        int cuty = cut;
        if (map.cols() == 1)
            cuty = 0;

        long naxes[2] =	{ long(map.rows() - 2 * cutx), long(map.cols() - 2 * cuty) };
                        


        // declare auto-pointer to FITS at function scope. Ensures no resources
        // leaked if something fails in dynamic allocation.
        auto_ptr<FITS> pFits(0);
        
        //remove file, otherwise get error from fits
        remove(fileName.c_str());
        pFits.reset(new FITS(fileName, DOUBLE_IMG, naxis, naxes));

        long nelements(1);
            
        nelements = naxes[0]*naxes[naxis];
        nelements = accumulate(&naxes[0], &naxes[naxis], 1, multiplies<long> ());

            
        valarray<float> array(nelements);

        for (int i = cutx; i < map.rows() - cutx; i++)
            for (int j = cuty; j < map.cols() - cuty; j++)
            {
                //to map it in the same way into the FITS as it is shown on the screen (readout is a the bottom)
                array[(i - cutx) + (j - cuty) * (map.rows() - 2 * cutx)] = (map(i, j)); 
            }

        pFits->pHDU().addKey("EXPOSURE", exposureTime, "Total Exposure Time");
        pFits->pHDU().write(1, nelements, array);

        pFits.reset();

}
//==============================================================================




//==============================================================================
/**
 * Wrapper class to handle fits files. Provides convenience methods to write fits files from a Blitz array of doubles.
 * 
 * @param exposureTime Exposure time, to be written in the FITS file header.
 * @param cut Controls how much of the edge of the image is not written to FITS file (should be half size of PSF).
 * @param &map Blitz Array<double, 2> to be written in the FITS file. 
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 */
void FileUtilities::FileUtilitiesWriteFITS(double exposureTime, int cut, Array<double, 2> &map, std::string fileName)
{
         
        LogManager::log << "    Creating FITS file " ;
        GlobalVariables::logManager.LogManagerShowLog();        

        
        fileName += ".fits";
        long naxis = 2;
        int cutx = cut;
        int cuty = cut;
        if (map.cols() == 1)
            cuty = 0;

        long naxes[2] =	{ long(map.rows() - 2 * cutx), long(map.cols() - 2 * cuty) };
                        


        // declare auto-pointer to FITS at function scope. Ensures no resources
        // leaked if something fails in dynamic allocation.
        auto_ptr<FITS> pFits(0);
        
        //remove file, otherwise get error from fits
        remove(fileName.c_str());
        pFits.reset(new FITS(fileName, DOUBLE_IMG, naxis, naxes));

        long nelements(1);
            
        nelements = naxes[0]*naxes[naxis];
        nelements = accumulate(&naxes[0], &naxes[naxis], 1, multiplies<long> ());

            
        valarray<float> array(nelements);

        for (int i = cutx; i < map.rows() - cutx; i++)
            for (int j = cuty; j < map.cols() - cuty; j++)
            {
                //to map it in the same way into the FITS as it is shown on the screen (readout is at the bottom)
                array[(i - cutx) + (j - cuty) * (map.rows() - 2 * cutx)] = (map(i, j)); 
            }

        pFits->pHDU().addKey("EXPOSURE", exposureTime, "Total Exposure Time");
        pFits->pHDU().write(1, nelements, array);

        pFits.reset();

}
//==============================================================================



//==============================================================================
/**
 * Wrapper class to handle fits files. Provides convenience methods to write fits files from a Blitz array of doubles.
 * 
 * @param exposureTime Exposure time, to be written in the FITS file header.
 * @param cut Controls how much of the edge of the image is not written to FITS file (should be half size of PSF).
 * @param numPrescanRows
 * @param numSmearingOverscanRows
 * @param biasRegisterMap Blitz Array<float, 2> containing the Smearing map.       
 * @param smearingMap Blitz Array<float, 2> containing the BIAS register map.  
 * @param &map Pointer to Blitz Array<double, 2> to be written in the FITS file. 
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 */
void FileUtilities::FileUtilitiesCutAndWriteFITS(double exposureTime, int cut, int numPrescanRows, int numSmearingOverscanRows, Array<float, 2>  biasRegisterMap, 
                                                 Array<float, 2>  smearingMap, Array<float, 2> &map, std::string fileName)
{
         
        LogManager::log << "    Creating FITS file " ;
        GlobalVariables::logManager.LogManagerShowLog();        

        
        fileName += ".fits";
        long naxis = 2;
        int cutx = cut;
        int cuty = cut;
        
         //Just to be sure that the map is not empty
        if (map.cols() == 0)
        {
            LogManager::log << "WARNING: This map is empty!";        
            GlobalVariables::logManager.LogManagerShowLog();         
        }
        //Just to be sure that the map is only a line and the cut is more than one line
        if (map.cols() == 1)
        {
            cuty = 0;
        }
        
        //Generate a new cutMap from the input map cropping the input value"cut" from up, down, right and left sides.
        Array<float, 2> cutMap(map.rows() - 2 * cutx, map.cols() - 2 * cuty);
        cutMap = map(Range(cutx, map.rows() - cutx - 1), Range(cuty, map.cols() - cuty - 1));

        //Generate a newMap from the cropped one but adding (numPrescanRows + numSmearingOverscanRows) rows
        Array<float, 2> newMap(cutMap.rows(), cutMap.cols() + numPrescanRows + numSmearingOverscanRows);
        
        //Add the BIAS register to the newMap from the first row to the numPrescan row
        newMap(Range::all(), Range(0, numPrescanRows - 1)) = biasRegisterMap(Range(cutx, map.rows() - 1), Range(0, numPrescanRows - 1));

        //Add the cropped subfield map to the newMap, starting in the numPrescan row.
        newMap(Range::all(), Range(numPrescanRows, numPrescanRows + cutMap.cols() - 1)) = cutMap;

        //Add the smearing register map to the newMap, starting at the end of the croppped subfield map.
        newMap(Range::all(), Range(numPrescanRows + cutMap.cols(), toEnd)) = smearingMap(Range(cutx, map.rows() - 1), Range(0, numSmearingOverscanRows - 1));
       
        
        
        
        long naxes[2] =	{ long(newMap.rows()), long(newMap.cols()) };
                        


        // declare auto-pointer to FITS at function scope. Ensures no resources
        // leaked if something fails in dynamic allocation.
        auto_ptr<FITS> pFits(0);
        
        //remove file, otherwise get error from fits
        remove(fileName.c_str());
        pFits.reset(new FITS(fileName, DOUBLE_IMG, naxis, naxes));

        long nelements(1);
            
        nelements = naxes[0]*naxes[naxis];
        nelements = accumulate(&naxes[0], &naxes[naxis], 1, multiplies<long> ());

            
        valarray<float> array(nelements);

        for (int i = 0; i < newMap.rows(); i++)
            for (int j = 0; j < newMap.cols(); j++)
            {
                //to map it in the same way into the FITS as it is shown on the screen (readout is at the bottom)
                array[(i) + (j) * (newMap.rows())] = int(round(newMap(i, j))); 
            }

        pFits->pHDU().addKey("EXPOSURE", exposureTime, "Total Exposure Time");
        pFits->pHDU().write(1, nelements, array);

        pFits.reset();
       
        LogManager::log << "    CREATED FITS file " ;
        GlobalVariables::logManager.LogManagerShowLog();        

}
//==============================================================================



//==============================================================================
/**
 * Wrapper class to read fits files. Provides convenience methods to read fits files.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param &flux Blitz Array<int, 2> to be read from the FITS file.
 */
void FileUtilities::FileUtilitiesReadFits( std::string fileName, Array<int, 2> &flux )
{

    if ( fileName == "" ) return;
    auto_ptr<FITS> pInfile ( new FITS ( fileName, Read, std::string("SCI") ) );
    PHDU& image = pInfile->pHDU();

    valarray<unsigned int>  contents;

    // read all user-specifed, coordinate, and checksum keys in the image
    image.readAllKeys();
    image.read ( contents );

    
    //dimensions of the image
    int ax1 ( image.axis ( 0 ) );
    int ax2 ( image.axis ( 1 ) );
      

    for (int k = 0; k < ax2; k+=10)
    {    
        ostream_iterator<short> c(std::cout,"\t");
        copy(&contents[k*ax1],&contents[(k+1)*ax1-1],c);
        std::cout <<&contents[k*ax1]<<"  "<<&contents[(k+1)] << '\n';
        
    }

    //completing the flux map with the FITS contents
    flux.resize ( ax1, ax2 );
    flux = 0.0;

    for ( int i = 0;i < ax1;i++ )
    {
        for ( int j = 0;j < ax2;j++ )
        {
            flux ( i, j ) = float ( contents[i+j*ax1] );
        }
    }

    pInfile.reset();
}
//==============================================================================



//==============================================================================
/**
 * Wrapper class to read fits files. Provides convenience methods to read fits files.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param &flux Blitz Array<double, 2> to be read from the FITS file.
 */
void FileUtilities::FileUtilitiesReadFits( std::string fileName, Array<double, 2> &flux )
{
    if ( fileName == "" ) return;
   
    auto_ptr<FITS> pInfile ( new FITS ( fileName, Read, true ) );
    PHDU& image = pInfile->pHDU();
    
    valarray<unsigned long>  contents;
    
    // read all user-specifed, coordinate, and checksum keys in the image
    image.readAllKeys();
    
    image.read ( contents );
    

    
    //dimensions of the image
    int ax1 ( image.axis ( 0 ) );
    int ax2 ( image.axis ( 1 ) );
  
    flux.resize ( ax1, ax2 );
    flux = 0.0;
    
    for ( int i = 0;i < ax1;i++ )
    {
        for ( int j = 0;j < ax2;j++ )
        {
            flux ( i, j ) = double ( contents[i+j*ax1] );
        }
    }
    
    pInfile.reset();
}
//==============================================================================



//==============================================================================
/**
 * Wrapper class to read fits files. Provides convenience methods to read fits files.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param &flux Blitz Array<float, 2> to be read from the FITS file.
 */
void FileUtilities::FileUtilitiesReadFits( std::string fileName, Array<float, 2> &flux )
{

    if ( fileName == "" )
    {
        std::cerr<<"\nError (FileUtilities::FileUtilitiesReadFits): Empty file name."<<std::endl;
        return;
    }
    
 

    std::auto_ptr<FITS> pInfile ( new FITS ( fileName, Read,  true ) );
       
    PHDU& image = pInfile->pHDU();
        

  
//    if(image.bitpix() != -32) 
//    {
//        std::cerr << "\nError (FileUtilities::FileUtilitiesReadFits): Fits image must be of type FLOAT. Make sure you are indicating the [SCI] extension in the input file name." << std::endl;
//        exit(1);
//    }
       

    valarray<unsigned long>  contents;

    // read all user-specifed, coordinate, and checksum keys in the image
    image.readAllKeys();
    image.read(contents);
       

    //dimensions of the image
    long nax(image.axes());
    long ax1(image.axis(0));
    long ax2(image.axis(1));
    long ax3(image.axis(2));


//    for (long k = 0; k < ax2; k+=10)
//    {    
//        ostream_iterator<short> c(std::cout,"\t");
//        copy(&contents[k*ax1],&contents[(k+1)*ax1-1],c);
//        //std::cout <<&contents[k*ax1]<<"  "<<&contents[(k+1)] << '\n';
//        
//    }

    
    //completing the flux map with the FITS contents
    flux.resize ( ax1, ax2 );
    flux = 0.0;
     
    for ( int i = 0;i < ax1;i++ )
    {
        for ( int j = 0;j < ax2;j++ )
        {
            flux ( i, j ) = float ( contents[i+j*ax1] );
        }
    }
    pInfile.reset();
        
}
//==============================================================================




//==============================================================================
/**
 * Wrapper class to read FITS files from an external source.
 * This will be used to get a raw image as input in order to add noise to that
 * image. This can be used for the WUVS Simulations, taking echelle spectra images
 * and include noise to them.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param &flux Blitz Array<float, 2> to be read from the FITS file.
 */
void FileUtilities::FileUtilitiesReadExternalFits( std::string fileName, Array<float, 2> &flux )
{
    
    FITS::setVerboseMode(false);

    // dimensions
    int ax1; // x
    int ax2; // y
    int ax3; // lambda
    int ax4; // #images

    int rv; // return value


    try
    {

        std::auto_ptr<FITS> pInfile(new FITS(fileName,Read,std::string("SCI")));
        ExtHDU& table = pInfile->extension("SCI");

        // read all the keywords, excluding those associated with columns.
        table.readAllKeys();

        int naxis=table.axes() ;

        ax1=table.axis(0) ; // x
        ax2=table.axis(1) ; // y


        if (naxis>=3) 
        {
            ax3=table.axis(2);
        }

        if (naxis>=4) 
        {
            ax4=table.axis(3); 
        }

        // reading image
        std::valarray<double>  image2D(ax1*ax2); // 2 dimension image
        table.read(image2D); // read image


        //completing the flux map with the FITS contents
        flux.resize ( ax1, ax2 );
        flux = 0.0;

        // Now, the next code alway work and does not have any drawbak
        int i=0;
        for (int r=0;r<ax1;r++)
        {      
          for (int c=0;c<ax2;c++)
            {
              flux ( r, c )=(float)image2D[i];
              i++;
            }
        }
        rv=0; // ok

    }

    catch (FitsException&) 
    // will catch all exceptions thrown by CCfits, including errors
    // found by cfitsio (status != 0)
    {

        std::cerr << " Fits Exception Thrown by FitsImage class \n";    
        std::cerr << " Fits file name : " << fileName << std::endl;
        rv=1; // problem   

    }
}
//==============================================================================





//==============================================================================
/**
 * Wrapper class to read fits files. Provides convenience methods to read fits files.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param map Blitz Array<double, 2> to be written to FITS file.
 */
void FileUtilities::FileUtilitiesWriteFits( std::string fileName, Array<double, 2> map, int cut )
{
    //qDebug ( "void MyFits::write ( std::string fileName, Array<double, 2> map, int cut )" );
    
    fileName += ".fits";
    long naxis =   2;
    int cutx = cut;
    int cuty = cut;
    if ( map.cols() == 1 ) cuty = 0;
    
    long naxes[2] = { long ( map.rows() - 2*cutx ), long ( map.cols() - 2*cuty ) };
    
    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    auto_ptr<FITS> pFits ( 0 );
    
    //remove file, otherwise get error from fits
    remove ( fileName.c_str() );
    pFits.reset ( new FITS ( fileName, ULONG_IMG , naxis , naxes ) );
    
    long nelements ( 1 );
    
    nelements = accumulate ( &naxes[0], &naxes[naxis], 1, multiplies<long>() );
    
    vector<long> extAx ( long ( map.rows() - 2*cutx ), long ( map.cols() - 2*cuty ) );
    
    valarray<int> array ( nelements );
    
    for ( int i = cutx;i < map.rows() - cutx;i++ )
        for ( int j = cuty;j < map.cols() - cuty;j++ )
            array[ ( i-cutx ) + ( j-cuty ) * ( map.rows()-2*cutx ) ] = long ( round ( map ( i, j ) ) ); //to map it in the same way into the FITS as I have it on the screen (readout is a the bottom )
    
    pFits->pHDU().write ( 1, nelements, array );
    
    pFits.reset();
    extAx.clear();
}
//==============================================================================




//==============================================================================
/**
 * Wrapper class to write fits files. Provides convenience methods to write fits files.
 *
 * @param fileName name (including the path, excluding the .fits) for the FITS file.
 * @param map Blitz Array<float, 2> to be written to FITS file.
 */
void FileUtilities::FileUtilitiesWriteFits( std::string fileName, Array<float, 2> map, int cut )
{
    
    fileName += ".fits";
    long naxis =   2;
    int cutx = cut;
    int cuty = cut;
    if ( map.cols() == 1 ) cuty = 0;

    long naxes[2] = { long ( map.rows() - 2*cutx ), long ( map.cols() - 2*cuty ) };

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    auto_ptr<FITS> pFits ( 0 );
    //remove file, otherwise get error from fits
    remove ( fileName.c_str() );
    pFits.reset ( new FITS ( fileName, ULONG_IMG , naxis , naxes ) );
    
    long nelements ( 1 );

    nelements = accumulate ( &naxes[0], &naxes[naxis], 1, multiplies<long>() );

    vector<long> extAx ( long ( map.rows() - 2*cutx ), long ( map.cols() - 2*cuty ) );
    
    valarray<int> array ( nelements );

    for ( int i = cutx;i < map.rows() - cutx;i++ )
    {
        for ( int j = cuty;j < map.cols() - cuty;j++ )
        {
            array[ ( i-cutx ) + ( j-cuty ) * ( map.rows()-2*cutx ) ] = long ( round ( map ( i, j ) ) ); //to map it in the same way into the FITS as I have it on the screen (readout is a the bottom )
        }
    }
    pFits->pHDU().write ( 1, nelements, array );

    pFits.reset();
    extAx.clear();
}
//==============================================================================



