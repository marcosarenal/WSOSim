///////////////////////////////////////////////////////////
//  ParamsFlatField.cpp
//  Implementation of the Class ParamsFlatField
//  Created on:      05-Dic-2012 16:59:59 PM
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



#include "ParamsFlatField.h"



//==============================================================================
/**
 * Constructor method
 */
ParamsFlatField::ParamsFlatField(){}
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
ParamsFlatField::~ParamsFlatField(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/**
 * Calculate and sets the FlatField map in the Dataset to be used in the corresponding
 * processing step (StepFlatField).
 * Compute the Flat Field at pixel or sub-pixel level and save them in a map that
 * is applied to the pixel or subPixelMap respectively.
 * The Flat Field at pixel or sub-pixel level are assumed to follow a spatial 1/f-response.
 */
void ParamsFlatField::paramsFlatFieldcreateFFMap(DataSet &m_DataSet)
{
    //Pointing to the DataSet
    p_DataSet = &m_DataSet;
    
    //Retrieving parameters from DataSet
    subPixelsPerPixel = p_DataSet->datasetGetsubPixelsPerPixel();
    subFieldSizeX = p_DataSet->datasetGetsubFieldSizeX();
    subFieldSizeY = p_DataSet->datasetGetsubFieldSizeY();
    
    //Fractional peak-to-peak amplitude of the pixel non-uniform sensitivity response.
    flatfieldPixelNoise = p_DataSet->datasetGetflatfieldPixelNoise();
    
    //Fractional white noise component of the sub-pixels.
    flatfieldWhiteNoise = p_DataSet->datasetGetflatfieldWhiteNoise();
    
    //Flatfield intra pixel-width in % at edge of pixel with 5\% lower sensitivity [\% of pixel size, rounded up]
    flatfieldIntraPixelWidth = p_DataSet->datasetGetflatfieldIntraPixelWidth();
    
    
    edge = ceil(subPixelsPerPixel * flatfieldIntraPixelWidth/ 100.); //flatfieldIntraPixelWidth was given as a %
    
    //Checks whether the flatfield map should be computed at pixel or sub-pixel level.
    //if sub-pixel white noise and the intra-pixel width are null,
    //the flatfield map should be computed at normal level.
    if (flatfieldWhiteNoise == 0 && flatfieldIntraPixelWidth == 0)
    {
        ParamsFlatField::paramsFlatFieldcreatePixelFFMap();
        
        //Setting the generated FlatField map in the Dataset
        m_DataSet.datasetSetFlatFieldMap(flatfieldMap);
        
    }
    
    //Else, there must be computed the subPixelFFMap.
    else
    {
                    
        intraPixelSensitivity.resize(edge);
        intraPixelSensitivity = 0.0;
        
        for (int i = 0; i < edge; i++)
        {
            intraPixelSensitivity(i) = 0.95; // The input parameter flatfieldIntraPixelWidth is defined as the
            //  flatfield intrapixel-width at edge of pixel with 5% lower sensitivity [% of pixel size, rounded up]
        }
        
        
        ParamsFlatField::paramsFlatFieldcreateSubPixelFFMap();
            
        //Setting the generated subFlatField map in the Dataset
        m_DataSet.datasetSetsubFlatFieldMap(subflatfieldMap);
        
            
        //Create a Flat Field map the same size as the subfield
        flatfieldMap.resize(subFieldSizeX, subFieldSizeY);

        //Rebin the subFlatFieldMap from subpixel level to pixel level to be used in photometry
        for (int i = 0; i < subFieldSizeX; i++)
        {        
            for (int j = 0; j < subFieldSizeY; j++)
            {        
                flatfieldMap(i, j) = sum(subflatfieldMap(Range(i * subPixelsPerPixel, (i + 1) * subPixelsPerPixel - 1), 
                                                         Range(j * subPixelsPerPixel, (j + 1) * subPixelsPerPixel - 1)));
            }
        } 
    
        //Setting the generated FlatField map in the Dataset
        m_DataSet.datasetSetFlatFieldMap(flatfieldMap);
    }
   
    
    
    LogManager::log<< "    FlatField map created ";
    GlobalVariables::logManager.LogManagerShowLog();
    
}
//==============================================================================




//==============================================================================
/**
 * Compute the Flat Field at sub-pixel level and save them in a map that can be applied to the subPixelMap.
 * The Flat Field at normal pixel level are assumed to follow a spatial 1/f-response. The white noise component is
 * applied to the sub-pixels. The sensitivity at the edge of each pixel is lower than at the center of a pixel.
 */
void ParamsFlatField::paramsFlatFieldcreateSubPixelFFMap()
{
    
    
    NormalUnit<double, ranlib::MersenneTwister, ranlib::independentState> normal;
    normal.seed(DataSet::seedRNG);

    //Create a Flat Field map the same size as the subfield
    subflatfieldMap.resize(subFieldSizeX * subPixelsPerPixel, subFieldSizeY * subPixelsPerPixel);
    subflatfieldMap = 0.0;

    int newj, newi;
    int modulo;


    //Take the higher power of two smaller than the subfield in order to contain variations at all spatial frequencies.
    int dimpow2 = 2;
    int m = max(subFieldSizeX, subFieldSizeY);
    
    do
    {
        dimpow2 *= 2;
    } while (dimpow2 <= m);

    
    Array<float, 2> flat;
    //flat is at subpixel level, flatfieldMap is at pixel level

    //create flatfield image that is bigger than data image (only section relevant for data is then taken)
    try
    {
        flat.resize(dimpow2, dimpow2);
        flat = 0.0;

    } 
    catch (int e)
    {
        cerr << "\nError (ParamsFlatField::paramsFlatFieldcreateSubPixelFFMap()): Cannot allocate enough memory for copying the sub-pixel matrix." << endl;
        exit(1);
    }
    
    
    
    //Goes through all the pixels in the flat field and gives random values on that pixel. 
    //Then gives the same value to that and (n-1)**2 pixels around that one.
    //Then gives the same value to that and (n-1)**3 pixels around that one.
    //...
    //... up to the flat size = dimpow2.
    //This provides variations at all spatial frequencies.    
    double w = 2;

    for (int z = dimpow2 / 2; z >= 2; z /= 2)
    {
        for (int x = 0; x < w; x++)
        {
            for (int y = 0; y < w; y++)
            {
                flat(Range(x * z, (x + 1) * z - 1), Range(y * z, (y + 1) * z - 1)) += sqrt(w) * normal.random();
            }
        }

    w = w * 2;
    }

    //normalize flat, set to between 0 and 1, and multiply with peak-to-peak noise amplitude
    flat -= min(flat);
    flat = flat / max(flat) - 0.5;
    flat *= flatfieldPixelNoise;
    
    for (int i = 0; i < subFieldSizeX; i++)
    {
        for (int j = 0; j < subFieldSizeY; j++)
        {
            for (int h = 0; h < subPixelsPerPixel; h++)
            {
                for (int c = 0; c < subPixelsPerPixel; c++)
		{
                    newi = (i * subPixelsPerPixel) + h;
                    newj = (j * subPixelsPerPixel) + c;

                    if (h < edge)
                    {
                        subflatfieldMap(newi, newj) = (1 + flat(i, j) + normal.random() * flatfieldWhiteNoise) * (intraPixelSensitivity(h));
                    }
                    else if (c < edge)
                    {
                        subflatfieldMap(newi, newj) = (1 + flat(i, j) + normal.random() * flatfieldWhiteNoise) * (intraPixelSensitivity(c));
                    }
                    else if (c >= subPixelsPerPixel - edge)
                    {
                        modulo = (c + edge) % subPixelsPerPixel;
                        subflatfieldMap(newi, newj) = (1 + flat(i, j) + normal.random() * flatfieldWhiteNoise) * (intraPixelSensitivity(modulo));
                    }
                    else if (h >= subPixelsPerPixel - edge)
                    {
                        modulo = (h + edge) % subPixelsPerPixel;
                        subflatfieldMap(newi, newj) = (1 + flat(i, j) + normal.random() * flatfieldWhiteNoise) * (intraPixelSensitivity(modulo));
                    }
                    else
                    {
                        subflatfieldMap(newi, newj) = (1 + flat(i, j) + normal.random() * flatfieldWhiteNoise);
                    }
                }
            }
        }
    }
    
    
    //Free memory
    flat.free();
    
    
}
//==============================================================================



//==============================================================================
/**
 * Compute the Flat Field at normal pixel level and save them in a map that can be applied to the pixelMap.
 * The Flat Field at normal pixel level are assumed to follow a spatial 1/f-response.
 */
void ParamsFlatField::paramsFlatFieldcreatePixelFFMap()
{
    
    NormalUnit<double, ranlib::MersenneTwister, ranlib::independentState> normal;
    normal.seed(DataSet::seedRNG);

    //Create a Flat Field map the same size as the subfield
    flatfieldMap.resize(subFieldSizeX, subFieldSizeY);

    //Take the higher power of two smaller than the subfield in order to contain variations at all spatial frequencies.
    int dimpow2 = 2;
    int m = max(subFieldSizeX, subFieldSizeY);
    
    do
    {
        dimpow2 *= 2;
    } while (dimpow2 <= m);

    //Create temporary flatfield image that is bigger than data image (only section relevant for data is then taken)
    Array<float, 2> flat;

    //Take a flat that have the pixel size of the higher power of two that is smaller than the subfield.
    flat.resize(dimpow2, dimpow2);
    flat = 0.;

    
    //Goes through all the pixels in the flat field and gives random values on that pixel. 
    //Then gives the same value to that and (n-1)**2 pixels around that one.
    //Then gives the same value to that and (n-1)**3 pixels around that one.
    //...
    //... up to the flat size = dimpow2.
    //This provides variations at all spatial frequencies.
    double w = 2;

    for (int z = dimpow2 / 2; z >= 2; z /= 2)
    {
        for (int x = 0; x < w; x++)
        {
            for (int y = 0; y < w; y++)
            {
                flat(Range(x * z, (x + 1) * z - 1), Range(y * z, (y + 1) * z - 1)) += sqrt(w) * normal.random();
            }
        }

        w = w * 2;
    }


    
    //Normalize flat, set to between 0 and 1, and multiply with peak-to-peak noise amplitude
    //Remove offset
    flat -= min(flat);
    
    //Divide every pixel by maximum value to normalize to 1.
    flat = flat / max(flat) - 0.5;
    flat *= (2 * flatfieldPixelNoise);

    //Copy the calculated flat field into the Flat Field Map, but only the subfield size
    flatfieldMap = 1. + flat(Range(0, subFieldSizeX - 1), Range(0, subFieldSizeY - 1));


    //Free memory
    flat.free();
    
}
//==============================================================================