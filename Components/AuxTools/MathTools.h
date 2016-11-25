///////////////////////////////////////////////////////////
//  MathTools.h
//  Math methods
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


#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <cmath>
#include <vector>
#include "blitz/array.h"

using namespace blitz;

class MathTools
{

public:
	MathTools();
	static void shakersort(Array<double, 1> &v);
	static void shakersort(Array<float, 1> &v);
	static void shakersort(vector<std::string> &v);
	static void shakersort(Array<double, 1> &v1, Array<double, 1> &v2);
	static void shakersort(Array<double, 1> &v1, Array<std::string, 1> &v2);
	static void shakersort(Array<double, 2> v, Array<double, 1> &vout, Array<int, 1> &index1, Array<int, 1> &index2);
	static void shakersort(Array<float, 2> v, Array<float, 1> &vout, Array<int, 1> &index1, Array<int, 1> &index2);
	static double median(Array<double, 1> v);
	static double median(Array<float, 1> v);
	static double interpolate(double x1, double y1, double x2, double y2, double x);
	static void interpolate(Array<double, 1> xmask, Array<double, 1> x, Array<double, 1> y, Array<double, 1> &temp);
	static int signof(double a)
	{
		return (a == 0) ? 0 : (a < 0 ? -1 : 1);
	}
	static double closestValue(double n, Array<double, 1> v);
	static void closestValue(double n, Array<double, 1> v, double &cv, int &index);

	static void getLOSAngle(double xe, double ze, double xt, double yt, double zt, double &siwi);
	static void karth2Polar(double x, double y, double z, double &radius, double &theta, double &fi);
	static void karth2Polar(float x, float y, float z, float &radius, float &theta, float &fi);
	static void polar2Karth(double theta, double fi, double radius, double &X, double &Y, double &Z);
	static void polar2Karth(float theta, float fi, float radius, float &X, float &Y, float &Z);
	static void polarTrans2Karth(double fitime, double e, double radius, double theta, double fi, double &X, double &Y, double &Z);
	static void correctSphAnglesOutOfLimit(double theta, double fi, double &thetaout, double &fiout);
	static void trans3D2D(double xe, double ze, double X, double Y, double Z, double &A, double &B);
	static void crossProduct(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double &Xn,
			double &Yn, double &Zn);
	static void crossProduct(double x0, double y0, double z0, double x1, double y1, double z1, double &Xn, double &Yn, double &Zn);
	static void crossProduct(float x0, float y0, float z0, float x1, float y1, float z1, float &Xn, float &Yn, float &Zn);
	static double absoluteValue(double x, double y, double z);
	static double absoluteValue(float x, float y, float z);
	static double absoluteValue(double x, double y);
	static void vectorDiff(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3);
	static double derivative(double funcvalue, double funcvalue_dx, double dx);
	static void crossProductSph(double r1, double t1, double f1, double r2, double t2, double f2, double r3, double t3, double f3,
			double &x, double &y, double &z);
	static double distance(double x1, double y1, double x2, double y2);
	static double getCircleSegmentArea(double r, double phi);
	static double getCircleOverlapArea(double x1, double y1, double r1, double x2, double y2, double r2);
	static double getAngleTriangle(double a, double b, double c);
	static double getAngularDistanceSphere(double fi1, double theta1, double fi2, double theta2);
	static void trans2DPolar(double X0, double Y0, double eTh, double eFi, double &sTh, double &sFi);
	static void normalProjectionPointOnPlane(double xp, double yp, double zp, double a, double b, double c, double &xd, double &yd,
			double &zd);
	static void normalProjectionPointOnPlane(float xp, float yp, float zp, float a, float b, float c, float &xd, float &yd, float &zd);
	static void rotation(double &xp, double &yp, double &zp, double ux, double uy, double uz, double theta);
	static void rotation(float &xp, float &yp, float &zp, float ux, float uy, float uz, float theta);
	static void normalize(double &x, double &y, double &z);
	static void normalize(float &x, float &y, float &z);
	static double distancePointLine(double xp, double yp, double zp, double ux, double uy, double uz);
	static void coordinateTransformation(double xp, double yp, double zp, double ix, double iy, double iz, double jx, double jy, double jz,
			double kx, double ky, double kz, double &xt, double &yt, double &zt);
	static void getGnomonicProjection(double ra, double dec, double ra0, double dec0, double gamma, double &x, double &y);
	static void getGnomonicProjection(float ra, float dec, float ra0, float dec0, float gamma, float &x, float &y);
	static void getInverseGnomonicProjection(double x, double y, double ra0, double dec0, double gamma, double &ra, double &dec);
	static void getInverseGnomonicProjection(float x, float y, float ra0, float dec0, float gamma, float &ra, float &dec);
	static void rotate(double xIn, double yIn, double x0, double y0, double angle, double &xOut, double &yOut);
	static void rotate(float xIn, float yIn, float x0, float y0, float angle, float &xOut, float &yOut);
	static void rotate(Array<double, 2> posIn, double x0, double y0, double yaw, double pitch, double roll, Array<double, 2> &posOut);
	static void rotate(Array<double, 2> posIn, double x0, double y0, double yaw, double roll, Array<double, 1> &posOut);
	static void rotateMap(Array<double, 2> mapIn, double x0, double y0, double angle, Array<double, 2> &mapOut);
	static void translate(double &x, double &y, double dx, double dy);
	static void translate(float &x, float &y, float dx, float dy);
	static double getAngle(double ux, double uy, double uz, double vx, double vy, double vz);
	static double getAngle(double ux, double uy, double vx, double vy);
	static double getInnerProduct(double ux, double uy, double uz, double vx, double vy, double vz);
	static double getInnerProduct(double ux, double uy, double vx, double vy);
};
#endif /* MATHTOOLS_H_ */
