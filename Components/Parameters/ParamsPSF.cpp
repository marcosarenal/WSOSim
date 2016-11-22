///////////////////////////////////////////////////////////
//  ParamsPSF.cpp
//  Implementation of the Class ParamsPSF
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



#include "ParamsPSF.h"





//==============================================================================
/**
 * Constructor method
 */
ParamsPSF::ParamsPSF(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsPSF::~ParamsPSF(){}
//==============================================================================


//==============================================================================
//Functions:
//==============================================================================
/**
 * This function adds the PSF effect to the image. Fist, the PSF is rotated, then the PSF
 * is convolved
 */
void ParamsPSF::paramsPSFCalculation(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    psfRotation = p_DataSet->datasetGetpsfRotation();
    psfRotationAngle = p_DataSet->datasetGetpsfRotationAngle();
    psfOrientation = p_DataSet->datasetGetpsfOrientation();
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels();
      
    raOpticalAxis = p_DataSet->datasetGetraOpticalAxis();
    decOpticalAxis = p_DataSet->datasetGetdecOpticalAxis();
    raCenterSubField = p_DataSet->datasetGetraCenterSubField();
    declCenterSubField = p_DataSet->datasetGetdeclCenterSubField();  
    
    xOpticalAxis = p_DataSet->datasetGetxOpticalAxis(); 
    yOpticalAxis = p_DataSet->datasetGetyOpticalAxis();
    
    subFieldZeroX = p_DataSet->datasetGetsubFieldZeroX();
    subFieldZeroY = p_DataSet->datasetGetsubFieldZeroY();
    
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    outputPath = p_DataSet->datasetGetOutputPath();
    prefix = p_DataSet->datasetGetPrefix();
    
    //Checking the PSF input parameters
    ParamsPSF::paramsPSFcheckPSFParams();        

        
    //Calculate the angular distance
    angDist = Constants::RAD2DEG * MathTools::getAngularDistanceSphere(raOpticalAxis, decOpticalAxis, 
                                                                       raCenterSubField * Constants::DEG2RAD, declCenterSubField * Constants::DEG2RAD);					

    LogManager::log << "    Angular distance of the center of the sub-field from the optical axis: " << angDist << " degrees";
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
    
    //Initialization of the rotation angle
    rotationAngle = 0; 
    
    //Rotating the PSF
    if (psfRotation == PSF_ROTATION_OPTICAL_AXIS)
    {

        //Compute xOpticalAxis and yOpticalAxis relative to the center of the sub-field
        xOpticalAxis -= (subFieldZeroX + subFieldSizeX / 2.);
        yOpticalAxis -= (subFieldZeroY + subFieldSizeY / 2.);

        //Compute orientation vector in x,y of pre-computed PSF
        psfX = cos(psfOrientation * Constants::DEG2RAD);
        psfY = sin(psfOrientation * Constants::DEG2RAD);
        rotationAngle = Constants::RAD2DEG * MathTools::getAngle(psfX, psfY, xOpticalAxis, yOpticalAxis);

        //The function above only returns values between 0 and 180 degrees. This is corrected here to 0 to 360
        if (psfX * yOpticalAxis - xOpticalAxis * psfY > 0)
        {
            rotationAngle = 360 - rotationAngle;
        }

        rotationAngle = 360 - rotationAngle;
        LogManager::log << "    PSF is rotated by " << rotationAngle << " degrees to point towards center of optical axis.";
        GlobalVariables::logManager.LogManagerAppendLogAndShow();
    }

    else if (psfRotation == PSF_ROTATION_ARBITRARY)
    {
        rotationAngle = psfRotationAngle;
        LogManager::log << "    PSF is rotated by " << rotationAngle << " degrees to point towards center of optical axis.";
        GlobalVariables::logManager.LogManagerAppendLogAndShow();
    }
    
    //Calculation of the PSF mask
    ParamsPSF::paramsPSFComputeMask(m_DataSet);

}
//==============================================================================



 
//==============================================================================
/**
 * Check whether the input parameters in the XML input file are correct to perform 
 * the PSF computations.
 */
void ParamsPSF::paramsPSFcheckPSFParams()
{
    // Check PSFRotation parameter	
    if (psfRotation != 0 && psfRotation != 1 && psfRotation != 2)
    {
        cerr << "\nError (ParamsPSF::paramsPSFcheckPSFParams()): / "
        "PSFRotation value must be: 0 = No Rotation / "
                "1 = Towards Optical Axis / "
                "2 = Arbitrary Rotation by PSFRotationAngle." << endl;
        exit(1);
    }
        
}
//==============================================================================



        
//==============================================================================
/**
 * Compute the PSF mask that will be used for image convolution. This method must be called 
 * before any image convolution is performed (but only once for a simulation if the PSF does not change).
 * Depending on the input parameters, this will compute a simple 
 * 2D-Gaussian PSF, or read the PSF from a file.
 * If the PSF is location dependent and a location-look-up table has been provided, 
 * this method will select the appropriate PSF from the table.
 * From a look-up table (psfLocationFile) this method selects the PSF that is closest to the given angular distance.
 */
void ParamsPSF::paramsPSFComputeMask(DataSet &m_DataSet)
{
    //Retrieving parameters from the DataSet            
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    psfLocationDependent = p_DataSet->datasetGetpsfLocationDependent();
    useGauss = p_DataSet->datasetGetuseGauss();
    
    //If useGauss is set to TRUE:
    if (useGauss)
    {
        //Retrieving parameter from the DataSet     
        psfGaussFWHM = p_DataSet->datasetGetpsfGaussFWHM();
        
        //Call the function to calculate a gaussian PSF shaped 
        ParamsPSF::paramsPSFsetGaussian();
    }

    //If psfLocationDependent is set to TRUE
    else if (psfLocationDependent)
    {
        //Retrieve the PSF location file name from the DataSet.
        psfLocationFile = p_DataSet->datasetGetpsfLocationFile();
        
        //When the PSF is given from a grid, this function chooses the 
        //closest PSF and sets that PSF as the psfFileName to be read in paramsPSFreadFromFile 
        ParamsPSF::paramsPSFgetPSFFileFromGrid(); 
    }

    // ELSE, read the PSF from the selected file.
    else
    {
        //Retrieve the PSF file name from the DataSet.
        psfFileName = p_DataSet->datasetGetpsfFileName();

        //Call to read the PSF parameters from the psfFileName file
        ParamsPSF::paramsPSFreadFromFile();
    } 
    
    	 
    //Finally, set the calculated PSF map in the DataSet
    m_DataSet.datasetSetPSFMap(mask);


}
//==============================================================================




//==============================================================================
/**
 * This function sets a Gaussian PSF.
 */
void ParamsPSF::paramsPSFsetGaussian()
{

    //Retrieving fromm the DataSet    
    psfGaussFWHM = p_DataSet->datasetGetpsfGaussFWHM();
 
    psfNumPixels = p_DataSet->datasetGetpsfNumPixels();
    psfCenterX = p_DataSet->datasetGetpsfCenterX();
    psfCenterY = p_DataSet->datasetGetpsfCenterY();
 

    if (psfGaussFWHM < 1. / subPixelsPerPixel)
    {
        cerr << "\nError (ParamsPSF::paramsPSFsetGaussian()): The FWHM of the Gaussian PSF (" << psfGaussFWHM << ") must not be smaller than the size of a sub-pixel (" << 
        1./ subPixelsPerPixel << ").\n Please increase the size of the Gaussian FWHM or increase the number of sub-pixels." << endl;
        exit(1);
    }

    double fwhm = psfGaussFWHM * subPixelsPerPixel;
    double fak = 1. / (fwhm * fwhm * Constants::Pi2);
    double fak2 = 2 * fwhm * fwhm;
    	
    Array<float, 2> gaussMask;
    gaussMask.resize(psfNumPixels, psfNumPixels);
    gaussMask=0.0;
		
    	
    firstIndex i;
    secondIndex j;
    
    //The Gaussian mask is generated being the center at (psfCenterX, psfCenterY) 
    gaussMask = fak * exp(-(pow(i - psfCenterX + 0.5, 2.) + pow(j - psfCenterY + 0.5, 2.)) / fak2);

    //The whole matrix is shifted (added 0.5) such that the center of the PSF will lie on the cross section of four pixels.
    //If the center of the input psf is for instance at the array coordinates (3,5), then its psfCenter is at (3.5, 5.5).
    //Below, the matrix is shifted by (-0.5, -0.5). Then, the center will be located at the intrinsic pixel coordinates (3,5) 


    //determine new mask size. the center of the PSF is at (xCenter,yCenter)
    //and the mask must have odd size and be quadratic
    //the center of the PSF will be shifted and additional rows/colums added that have zero-flux
    //new size of PSF mask
    double xCenter = int(psfCenterX);
    double yCenter = int(psfCenterY);
    int newSize = max(max(xCenter + 1, yCenter + 1), max(psfNumPixels - xCenter, psfNumPixels - yCenter)) * 2 - 1;

    //generate a new mask with flux 0 and put the fileMask into it such that the center (x,y) is at the center ((xsize-1)/2) of the new mask
    Array<float, 2> newGaussMask(newSize, newSize);
    newGaussMask = 0.;
    double center = (newSize - 1) / 2;

    newGaussMask(Range(center - xCenter, center - xCenter + psfNumPixels - 1), Range(center - yCenter, center - yCenter + psfNumPixels - 1))
                = gaussMask;



    double xbinning = double(psfSubPixels) / double(subPixelsPerPixel);
    double ybinning = xbinning;

    int size = int(max(newSize / xbinning, newSize / ybinning));

    //rotate the PSF       
    if (psfRotation != NO_PSF_ROTATION)
    {
//        ParamsPSF::paramsPSFrotatePSF(newGaussMask, orientationAngle);
    }
        
    //rebin to subpixel level
    if (xbinning >= 1 && ybinning >= 1)
    {
        if (size % 2 == 0)
        {
            size++;
        }
        mask.resize(size, size);
        mask = 0.;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                {
                    mask(i, j) = sum(newGaussMask(Range(i * xbinning, (i + 1) * xbinning - 1), Range(j * ybinning, (j + 1) * ybinning - 1)));
                }
            }
        }
    }

    //Write the generated Gaussian PSF mask to ASCII file
    ParamsPSF::paramsPSFwriteSubPixelsToASCIIFile(outputPath+"/"+prefix+"/psfGaussianMask.dat");
 
    /** 
    //Patch to generate an output of the Gaussian PSFmap at pixel level.
    //This patch rebins the PSFmap to pixel level and writes it in a FITS file in order to assess
    //the Gaussian PSF generation function. 
    //IMPORTANT NOTE: It must be noted that FITS writing may create problems at this stage of the processing.
    //This piece of code should only be uncommented to generate the PSFmap at pixel level. 
    //The post-processing afterwards may generate errors.
    double xBins = double(size) / double(subPixelsPerPixel);
	double yBins = xBins;
	
    //generate a new newRebinnedPsf with flux 0 and put the fileMask into it such that the center (x,y) is at the center ((xsize-1)/2) of the new newRebinnedPsf
    Array<float, 2> newRebinnedPsf(xBins, xBins);

    int psfRebinnedSize = int(double(size) / double(xBins));


    newRebinnedPsf = 0.;
    
    for (int i = 0; i < xBins; i++)
    {
        for (int j = 0; j < yBins; j++)
        {
            if ((i + 1) * subPixelsPerPixel - 1 < size && (j + 1) * subPixelsPerPixel - 1 < size)
            {
                newRebinnedPsf(i, j) = sum(mask(Range(i * subPixelsPerPixel, (i + 1) * subPixelsPerPixel - 1), 
                                                Range(j * subPixelsPerPixel, (j + 1) * subPixelsPerPixel - 1)));
            
            }
        }
    }

  
    
    
    
    string outputPath = p_DataSet->datasetGetOutputPath();
    string    prefix = p_DataSet->datasetGetPrefix();
    string    outputDir = outputPath + "/" + prefix;    


        
    FileUtilities::FileUtilitiesWriteFITS(0, 0, newRebinnedPsf, outputDir + "/" + prefix + "newRebinnedPsf");

    newRebinnedPsf.free();
    */

    
    if (xbinning == 0.5 && ybinning == 0.5)
    {
            //mask must have odd size
            int size2 = newSize * 2 + 1;

            mask.resize(size2, size2);

            mask = 0;
            xbinning = ybinning = 1;

            for (int i = 0; i < size; i++)
    {
                    for (int j = 0; j < size; j++)
                    {
                            if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                            {
                                    mask(i * 2, j * 2) = sum(newGaussMask(Range(i * xbinning, (i + 1) * xbinning - 1), Range(j * ybinning, (j + 1) * ybinning
                                                    - 1)));
                            }
                    }
    }
            for (int i = 0; i < size - 1; i++)
    {
                    for (int j = 0; j < size - 1; j++)
                    {
                            if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                            {
                                    mask(i * 2 + 1, j * 2) = (mask(i * 2, j * 2) + mask(i * 2 + 2, j * 2)) / 2.;
                                    mask(i * 2, j * 2 + 1) = (mask(i * 2, j * 2) + mask(i * 2, j * 2 + 2)) / 2.;
                                    mask(i * 2 + 1, j * 2 + 1) = (mask(i * 2, j * 2) + mask(i * 2 + 2, j * 2) + mask(i * 2, j * 2 + 2) + mask(i * 2 + 2, j
                                                    * 2 + 2)) / 4.;
                            }
                    }
    }

    }

        
    //Free masks
    gaussMask.free();
    newGaussMask.free();


}
//==============================================================================




//==============================================================================
/**
 * This function reads the PSF from the selected ASCII file.
 */
void ParamsPSF::paramsPSFreadFromFile()
{
    
    //Retrieving parameters from DataSet
    psfSubPixels = p_DataSet->datasetGetpsfSubPixels();
    psfNumPixels = p_DataSet->datasetGetpsfNumPixels();
        
    psfCenterX = p_DataSet->datasetGetpsfCenterX();
    psfCenterY = p_DataSet->datasetGetpsfCenterY();
    
    
    
    if (floor(double(psfNumPixels) / double(psfSubPixels)) != ceil(double(psfNumPixels) / double(psfSubPixels)))
    {
        std::cerr << "\nError (ParamsPSF::paramsPSFreadFromFile()): Size of PSF (" << psfNumPixels << "=Number of Rows) must be "
                        " an integer factor of the number of sub-pixels ("<< psfSubPixels << ")." << std::endl;
        cerr << "\nModify the parameters NumPixels and SubPixels of PSFParameters accordingly." << endl;
        exit(1);
    }

    if ((double(psfSubPixels) / double(subPixelsPerPixel) != round(double(psfSubPixels) / double(subPixelsPerPixel)) && 
        double(subPixelsPerPixel) / double(psfSubPixels) != 2) || 
        double(psfSubPixels) / double(subPixelsPerPixel) <2 )

    {
        cerr << "\nError (ParamsPSF::paramsPSFreadFromFile()): The number of sub-pixels (" << subPixelsPerPixel<< ")"
                          " must be a integer multiple of the number of PSF sub-pixels per pixel (" << psfSubPixels << ") ";
        cerr << "of the input PSF used in the simulation." << endl;
        cerr << "\nModify the SubPixels parameter in the input parameters file accordingly." << endl;
        exit(1);
    }

    Array<float, 2> fileMask;
    fileMask.resize(psfNumPixels, psfNumPixels);
    fileMask=0.0;

    //Read PSF file into array fileMask
    double dummy;
    string str, line;
    int lineNumber = 0;
    ifstream myfile(psfFileName.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            istringstream iss(line);

            for (int i = 0; i < psfNumPixels; i++)
            {
                iss >> dummy;
                fileMask(i, lineNumber) = dummy;
            }
            lineNumber++;

        }
        myfile.close();
    }
    else
    {
        cerr << "\nError (ParamsPSF::paramsPSFreadFromFile()):: Unable to open file " << psfFileName << endl;
        exit(1);
    }

    LogManager::log << "    Integrated intensity of input PSF : " << sum(fileMask);
    GlobalVariables::logManager.LogManagerAppendLogAndShow();
        

                
    //determine offset of the PSF center from the cross section with x < xcenter and y < ycenter
    double dx = psfCenterX - int(psfCenterX);
    double dy = psfCenterY - int(psfCenterY);

    //shift the whole matrix such that the center of the PSF will lie on the cross section of four pixels.
    //do this with bilinear interpolation
    //If the center of the input psf is for instance at the array coordinates (3,5), then its psfCenter is at (3.5, 5.5).
    //Below, the matrix is shifted by (-0.5, -0.5). Then, the center will be located at the intrinsic pixel coordinates (3,5) 
    double p1, p2, p3, p4;

    for (int i = 0; i < psfNumPixels; i++)
    {
        for (int j = 0; j < psfNumPixels; j++)
        {
            p1 = fileMask(i, j);
            if (i < psfNumPixels - 1)
                    p2 = fileMask(i + 1, j);
            else
                    p2 = 0.;
            if (j < psfNumPixels - 1)
                    p3 = fileMask(i, j + 1);
            else
                    p3 = 0.;
            if (i < psfNumPixels - 1 && j < psfNumPixels - 1)
                    p4 = fileMask(i + 1, j + 1);
            else
                    p4 = 0.;
            //bilinear interpolation
            fileMask(i, j) = p1 * (1 - dx) * (1 - dy) + p2 * dx * (1 - dy) + p3 * (1 - dx) * dy + p4 * dx * dy;
        }
    }

    //determine new mask size. the center of the PSF is at (xCenter,yCenter)
    //and the mask must have odd size and be quadratic
    //the center of the PSF will be shifted and additional rows/colums added that have zero-flux
    //new size of PSF mask
    double xCenter = int(psfCenterX);
    double yCenter = int(psfCenterY);
    int newSize = max(max(xCenter + 1, yCenter + 1), max(psfNumPixels - xCenter, psfNumPixels - yCenter)) * 2 - 1;
    //generate a new mask with flux 0 and put the fileMask into it such that the center (x,y) is at the center ((xsize-1)/2) of the new mask
    Array<float, 2> newMask(newSize, newSize);
    newMask = 0.;
    double center = (newSize - 1) / 2;

    newMask(Range(center - xCenter, center - xCenter + psfNumPixels - 1), Range(center - yCenter, center - yCenter + psfNumPixels - 1))
                    = fileMask;

    double xbinning = double(psfSubPixels) / double(subPixelsPerPixel);
    double ybinning = xbinning;

    int size = int(max(newSize / xbinning, newSize / ybinning));

    //rotate the PSF       
    if (psfRotation != NO_PSF_ROTATION)
    {
//        ParamsPSF::paramsPSFrotatePSF(newMask, orientationAngle);
    }
        
    //rebin
    if (xbinning >= 1 && ybinning >= 1)
    {
        if (size % 2 == 0)
        {
            size++;
        }
        mask.resize(size, size);
        mask = 0.;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                {
                    mask(i, j) = sum(newMask(Range(i * xbinning, (i + 1) * xbinning - 1), Range(j * ybinning, (j + 1) * ybinning - 1)));
                }
            }
        }
    }

    ParamsPSF::paramsPSFwriteSubPixelsToASCIIFile(outputPath+"/"+prefix+"/psfReadMask.dat");


    if (xbinning == 0.5 && ybinning == 0.5)
    {
            //mask must have odd size
            int size2 = newSize * 2 + 1;

            mask.resize(size2, size2);

            mask = 0;
            xbinning = ybinning = 1;

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                    {
                            mask(i * 2, j * 2) = sum(newMask(Range(i * xbinning, (i + 1) * xbinning - 1), Range(j * ybinning, (j + 1) * ybinning
                                            - 1)));
                    }
                }
            }
            for (int i = 0; i < size - 1; i++)
            {
                for (int j = 0; j < size - 1; j++)
                {
                    if ((i + 1) * xbinning - 1 < newSize && (j + 1) * ybinning - 1 < newSize)
                    {
                        mask(i * 2 + 1, j * 2) = (mask(i * 2, j * 2) + mask(i * 2 + 2, j * 2)) / 2.;
                        mask(i * 2, j * 2 + 1) = (mask(i * 2, j * 2) + mask(i * 2, j * 2 + 2)) / 2.;
                        mask(i * 2 + 1, j * 2 + 1) = (mask(i * 2, j * 2) + mask(i * 2 + 2, j * 2) + mask(i * 2, j * 2 + 2) + mask(i * 2 + 2, j
                                        * 2 + 2)) / 4.;
                    }
                }
            }

    }

    LogManager::log << "    Integrated intensity of processed PSF : " << sum(mask);
    GlobalVariables::logManager.LogManagerAppendLogAndShow();

    //Free masks
    fileMask.free();
    newMask.free();
    
}
//==============================================================================




//==============================================================================
void ParamsPSF::paramsPSFgetPSFFileFromGrid()
{
    
    LogManager::log << "Reading the file with PSF locations ... " << psfLocationFile << endl;
    GlobalVariables::logManager.LogManagerShowLog(); 
    
	ifstream myfile(psfLocationFile.c_str());
	if (!myfile.is_open())
	{
		cerr << "\nError (ParamsPSF::paramsPSFgetPSFFileFromGrid()): Unable to open file " << psfLocationFile << endl;
		exit(1);
	}

	//count the number of lines in the psfLocationFile to know how big the array has to be
	int numLines = FileUtilities::countLines(psfLocationFile);
	if (numLines == 0)
	{
		cerr << "\nError (ParamsPSF::paramsPSFgetPSFFileFromGrid()): File " << psfLocationFile << " is empty." << endl;
		exit(1);
	}

	Array<double, 1> fdist(numLines);
	Array<string, 1> fname(numLines);

	string n;
	double d;
	int i = 0;

	//read the file names and distances from the psfLocationFile into the arrays
	while (myfile >> n >> d)
	{
            fname(i) = n;
            if (!FileUtilities::fileExists(fname(i)))
            {
                    cerr << "\nError (ParamsPSF::paramsPSFgetPSFFileFromGrid()): Unable to open file " << fname(i) << endl;
                    exit(1);
            }
            fdist(i++) = d;

	}
	myfile.close();

	//determine, which psf is closest to the input angular distance
	//first sort both arrays according to their angular distance
	MathTools::shakersort(fdist, fname);

	double closestValue;
	int index;
	MathTools::closestValue(angDist, fdist, closestValue, index);
	psfFileName = fname(index);

	GlobalVariables::logManager.log
			<< "Pre-computed PSF whose angular distance is closest to the center of the sub-field has following parameters:";
	GlobalVariables::logManager.LogManagerAppendLog();
	LogManager::log<< "File name: " << psfFileName << "(Distance to optical axis: " << closestValue << ")";
	GlobalVariables::logManager.LogManagerAppendLog();

}
//==============================================================================




//==============================================================================
/**
 * This function writes the sub-pixels PSF mask to an ASCII file. This PSF mask will be used forehand in 
 * the processing.
 * @param fileName
 */
void ParamsPSF::paramsPSFwriteSubPixelsToASCIIFile(string fileName)
{

    //Open the PSF file
	ofstream out(fileName.c_str());
	if (!out.is_open())
	{
		cerr << "\nError (ParamsPSF::paramsPSFwriteSubPixelsToASCIIFile()):Unable to open output PSF file";
		exit(1);
	}
	out.precision(10);
     

    //Write PSF to file
	for (int i = 0; i < mask.rows(); i++)
	{
		for (int j = 0; j < mask.cols(); j++)
        {
			out << (i - (mask.rows() - 1.) / 2.) / double(subPixelsPerPixel) << " " << (j - (mask.cols() - 1.) / 2.)
					/ double(subPixelsPerPixel) << " " << mask(i, j) << " " << endl;
        }
		out << endl;
	}
}
//==============================================================================




//==============================================================================
/**
 * Rotate the PSF mask around the point (psfCenterX,psfCenterY) using the library OpenCV.
 * This algorithm has been originally implemented by J.J. Green and adapted for PLATOSim
 * First convert the file mask to the opencv format
 * @param mask
 * @param orientationAngle
 */
void ParamsPSF::paramsPSFrotatePSF(Array<float, 2> &mask, double orientationAngle)
{
//
//	int n = mask.rows();
//
//	CvSize size = { n, n };
//	IplImage *src = cvCreateImage(size, IPL_DEPTH_32F, 1);
//
//	for (int i = 0; i < n; i++)
//    {
//		for (int j = 0; j < n; j++)
//        {
//			cvSetReal2D(src, i, j, mask(i, j));
//        }
//    }
//
//	//Perform the rotation with opencv
//	CvMat* Mrot = cvCreateMat(2, 3, CV_32FC1);
//	IplImage *dst;
//
//	CvPoint2D32f centre = cvPoint2D32f(n / 2. - 1, n / 2. - 1);
//
//	cv2DRotationMatrix(centre, orientationAngle, 1.0, Mrot);
//
//	dst = cvCloneImage(src);
//	dst->origin = src->origin;
//	cvZero(dst);
//
//	CvScalar fill = cvScalarAll(0.0);
//
//	cvWarpAffine(src, dst, Mrot, CV_INTER_CUBIC, fill);
//
//	cvReleaseImage(&src);
//	cvReleaseMat(&Mrot);
//
//	//Convert the rotated mask back to a blitz-Array
//	for (int i = 0; i < size.width; i++)
//    {
//		for (int j = 0; j < size.height; j++)
//        {
//			mask(i, j) = cvGetReal2D(dst, i, j);
//        }
//    }
//
//	cvReleaseImage(&dst);

}
//==============================================================================
