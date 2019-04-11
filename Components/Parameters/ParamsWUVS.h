///////////////////////////////////////////////////////////
//  ParamsWUVS.h
//  Implementation of the Class ParamsWUVS
//  Created on:      12-June-2015 1:59:58 PM
//  Original author: pablo
///////////////////////////////////////////////////////////
/*
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



#ifndef ParamsWUVS_H_
#define ParamsWUVS_H_

#include "blitz/array.h"
#include "DataSet.h"
#include "Constants.h"
#include "MathTools.h"


using namespace blitz;


/**
 * Defines all properties of the CCD, contains the sub- and normal-pixel maps and
 * contains several methods and algorithms for dealing with a CCD. Describes the
 * physical and electrical properties of the CCD. In the simulations, it is
 * assumed that the CCDs of all telescopes have the same general properties. The
 * read-out direction of the CCD is assumed to be oriented in negative y-direction
 * and the read-out strip is below the y=0 row. The field of view of the CCD
 * determines which stars affect the sub-field through read-out smearing (for a
 * charge transfer time > 0 s). The flatfield ( pixel non-uniform response) is
 * computed by considering a spatial 1/f-response of the sensitivity. Due to the
 * consideration of computing time, not the complete CCD is modeled with sub-
 * pixel precision but only a small sub-field with a side length of a few hundred
 * pixels.
 */
class ParamsWUVS
{

public:

	ParamsWUVS();
        virtual ~ParamsWUVS();
        void     ParamsWUVScalculation(DataSet &m_DataSet);
        void     ParamsWUVSresizeCCD(DataSet &m_DataSet);
        void     ParamsWUVScheckCCDParams(DataSet &m_DataSet);
        void     ParamsWUVSmappingCCDParams();
        void     ParamsWUVSsubmappingCCDParams();
        void     ParamsWUVSPixelOffsetToRad();
        void     ParamsWUVScomputeCCDOrientation();
        void     ParamsWUVScomputeCCDOrientation(double raJitterAxis, double declJitterAxis, 
                                                double yaw, double pitch, double roll);
        void     ParamsWUVScomputeRadiusFOVandCenter(); 
        void     ParamsWUVSgetTransformationCCDToSky(double xStarCCD, double yStarCCD, 
                                                    double &raS, double &decS);


private:
        
        DataSet*    p_DataSet;                         //Pointer to the DataSet to retrieve parameters from it.                                        

        int     ccdSizeX, ccdSizeY, subFieldZeroX, subFieldZeroY;       //Parameter retrieved from DataSet.
        int     subFieldSizeX, subFieldSizeY               ;            //Parameter retrieved from DataSet.
        double  pixelSize;                                              //Parameter retrieved from DataSet.
        double  originOffsetXmm, originOffsetYmm, pixelScale;           //Parameter retrieved from DataSet.
        double  ccdOriginOffsetX, ccdOriginOffsetY;                     //Parameter retrieved from DataSet.
        double  ccdOrientation;                                         //Parameter retrieved from DataSet.
        double  areaTelescope, transEff, quantEff, gain;                //Parameter retrieved from DataSet.
        int     fullWellSat, digitalSat;                                //Parameter retrieved from DataSet.
        
        int     edgePixels;                                             //Parameter retrieved from DataSet.
        int     halfPSFSize;
        int     psfSize, subPixelsPerPixel, psfSubPixels;               //Parameter retrieved from DataSet.
        double  opticalAxisRACenter, opticalAxisDecCenter, orientationFocalPlane;//Parameter retrieved from DataSet.
        double  xFOVSubField, yFOVSubField, radiusFOVCCD, raCenterSubField, declCenterSubField;//Parameter retrieved from DataSet.

        double  raOpticalAxis, decOpticalAxis, rotationAngleOA;         //Updated parameter to be set into the DataSet.
        double  meanCTE, fluxm0, readOutNoise, readOutTime;      //Parameter retrieved from DataSet.
        double  electronicOffset, exposureTime;                      //Parameter retrieved from DataSet.
        int     numPrescanRows, numSmearingOverscanRows;                //Parameter retrieved from DataSet.

        double flatfieldWhiteNoise, flatfieldPixelNoise, cosmicHitRate; //Parameter retrieved from DataSet.
        
        Array<float, 2> initPixelMap;                             //Blitz map set into the DataSet.
        Array<float, 2> initSubPixelMap;                          //Blitz map set into the DataSet.        
        Array<float, 2> pixelMap;                                 //Blitz map set into the DataSet.
        Array<float, 2> subPixelMap;                              //Blitz map set into the DataSet.


};                                                                      //Map array retrieved from DataSet.
#endif /* ParamsWUVS_H_ */
