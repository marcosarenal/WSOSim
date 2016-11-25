///////////////////////////////////////////////////////////
//  MathTools.cpp
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


#include "MathTools.h"
#include "Constants.h"
#include "Statistics.h"



//==============================================================================
/**
 * Constructor method
 */
MathTools::MathTools(){}
//==============================================================================




//==============================================================================
//Functions:
//==============================================================================
/*
//==============================================================================
 void MathTools::sort ( QVector<int> &v1, QVector<int> &v2 )
 {
 if ( v1.size() != v2.size() ) return;

 int help1;
 int help2;
 for ( int i = v1.size() - 1;i > 0;i-- )
 {
 for ( int j = 0;j < i;	j++ )
 {
 if ( v1[j] > v1[j+1] )
 {
 help1 = v1[j];
 v1[j] = v1[j+1];
 v1[j+1] = help1;

 help2 = v2[j];
 v2[j] = v2[j+1];
 v2[j+1] = help2;
 }
 }
 }
 }
//==============================================================================




//==============================================================================
 void MathTools::sort ( QVector<double> v1, QVector<double> &v2 )
 {
 if ( v1.size() != v2.size() ) return;

 double help1, help2;
 for ( int i = v1.size() - 1;i > 0;i-- )
 {
 for ( int j = 0;j < i;	j++ )
 {
 if ( v1[j] > v1[j+1] )
 {
 help1 = v1[j];
 v1[j] = v1[j+1];
 v1[j+1] = help1;

 help2 = v2[j];
 v2[j] = v2[j+1];
 v2[j+1] = help2;
 }
 }
 }
 }
//==============================================================================




//==============================================================================
 void MathTools::sort ( QVector<double> v1, QList<QString> &v2 )
 {
 if ( v1.size() != v2.size() ) return;

 double help1;
 QString help2;
 for ( int i = v1.size() - 1;i > 0;i-- )
 {
 for ( int j = 0;j < i;	j++ )
 {
 if ( v1[j] > v1[j+1] )
 {
 help1 = v1[j];
 v1[j] = v1[j+1];
 v1[j+1] = help1;

 help2 = v2[j];
 v2[j] = v2[j+1];
 v2[j+1] = help2;
 }
 }
 }
 }

//==============================================================================




//==============================================================================
 QVector<double> MathTools::sort ( QVector<double> v1 )
 {
 QVector<double> v = v1;
 double help1;

 for ( int i = v.size() - 1;i > 0;i-- )
 {
 for ( int j = 0;j < i;	j++ )
 {
 if ( v[j] > v[j+1] )
 {
 help1 = v[j];
 v[j] = v[j+1];
 v[j+1] = help1;
 }
 }
 }

 return v;
 }



//==============================================================================




//==============================================================================
 void MathTools::sorting3D ( QVector<double> &x1, QStringList &x2, QVector<double> &x3 )
 {
 double hilf1, hilf3;
 QString hilf2;

 for ( int i = x1.size() - 1;i > 0;i-- )
 {
 for ( int j = 0;j < i;j++ )
 {
 if ( x1[j] > x1[j+1] )
 {
 hilf1 = x1[j];
 x1[j] = x1[j+1];
 x1[j+1] = hilf1;

 hilf2 = x2[j];
 x2[j] = x2[j+1];
 x2[j+1] = hilf2;

 hilf3 = x3[j];
 x3[j] = x3[j+1];
 x3[j+1] = hilf3;
 }
 }
 }
 }

 */
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 */
void MathTools::shakersort(Array<double, 1> &v)
{
	int n = v.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (v(i - 1) > v(i))
			{
				swap(v(i - 1), v(i));
				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (v(i - 1) > v(i))
			{
				swap(v(i - 1), v(i));
				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 */
void MathTools::shakersort(Array<float, 1> &v)
{
	int n = v.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (v(i - 1) > v(i))
			{
				swap(v(i - 1), v(i));
				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (v(i - 1) > v(i))
			{
				swap(v(i - 1), v(i));
				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 */
void MathTools::shakersort(vector<std::string> &v)
{
	int n = v.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (v[i - 1] > v[i])
			{
				swap(v[i - 1], v[i]);
				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (v[i - 1] > v[i])
			{
				swap(v[i - 1], v[i]);
				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);
}
//==============================================================================




//==============================================================================
/**
 * Sort the arrays v1 and v2 according to the increasing order of v1.
 * @param v1
 * @param v2
 */
void MathTools::shakersort(Array<double, 1> &v1, Array<double, 1> &v2)
{
	
	int n = v1.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (v1(i - 1) > v1(i))
			{
				swap(v1(i - 1), v1(i));
				swap(v2(i - 1), v2(i));

				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (v1(i - 1) > v1(i))
			{
				swap(v1(i - 1), v1(i));
				swap(v2(i - 1), v2(i));

				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);
}

//==============================================================================




//==============================================================================
/**
 * sort the arrays v1 and v2 according to the increasing order of v1
 * @param v1
 * @param v2
 */
void MathTools::shakersort(Array<double, 1> &v1, Array<std::string, 1> &v2)
{
	//
	int n = v1.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (v1(i - 1) > v1(i))
			{
				swap(v1(i - 1), v1(i));
				swap(v2(i - 1), v2(i));

				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (v1(i - 1) > v1(i))
			{
				swap(v1(i - 1), v1(i));
				swap(v2(i - 1), v2(i));

				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);
}

/*
//==============================================================================




//==============================================================================
 void MathTools::shakersort ( Array<double, 1> v, Array<double, 1> &vout )
 {
 vout.resize ( v.size() );
 int n = v.size();
 int lower = 1;
 int higher = n - 1;
 int lastswap = n - 1;

 // down up down up down up down up down up down ...
 do {

 // bubble down
 for ( int i = higher; i >= lower; i-- )
 {
 if ( v ( i - 1 ) > v ( i ) )
 {
 swap ( v ( i - 1 ), v ( i ) );
 lastswap = i;
 }
 }
 lower = lastswap + 1;

 // bubble up
 for ( int i = lower; i <= higher; i++ )
 {
 if ( v ( i - 1 ) > v ( i ) )
 {
 swap ( v ( i - 1 ), v ( i ) );
 lastswap = i;
 }
 }
 higher = lastswap - 1;

 }
 while ( lower <= higher );

 vout = v;
 }
 */
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 * @param vout
 * @param index1
 * @param index2
 */
void MathTools::shakersort(Array<double, 2> v, Array<double, 1> &vout, Array<int, 1> &index1, Array<int, 1> &index2)
{
	vout.resize(v.rows() * v.cols());
    vout = 0.0;
    
	index1.resize(vout.size());
	index1=0;
    
	index2.resize(vout.size());
    index2=0;

	for (int i = 0; i < v.rows(); i++)
		for (int j = 0; j < v.cols(); j++)
		{
			vout(i * v.cols() + j) = v(i, j);
			index1(i * v.cols() + j) = i;
			index2(i * v.cols() + j) = j;
		}

	int n = vout.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (vout(i - 1) > vout(i))
			{
				swap(vout(i - 1), vout(i));
				swap(index1(i - 1), index1(i));
				swap(index2(i - 1), index2(i));
				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (vout(i - 1) > vout(i))
			{
				swap(vout(i - 1), vout(i));
				swap(index1(i - 1), index1(i));
				swap(index2(i - 1), index2(i));
				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 * @param vout
 * @param index1
 * @param index2
 */
void MathTools::shakersort(Array<float, 2> v, Array<float, 1> &vout, Array<int, 1> &index1, Array<int, 1> &index2)
{
	vout.resize(v.rows() * v.cols());
    vout = 0.0;
    
	index1.resize(vout.size());
	index1 = 0;
    
	index2.resize(vout.size());
    index2 = 0;

	for (int i = 0; i < v.rows(); i++)
		for (int j = 0; j < v.cols(); j++)
		{
			vout(i * v.cols() + j) = v(i, j);
			index1(i * v.cols() + j) = i;
			index2(i * v.cols() + j) = j;
		}

	int n = vout.size();
	int lower = 1;
	int higher = n - 1;
	int lastswap = n - 1;

	// down up down up down up down up down up down ...
	do
	{

		// bubble down
		for (int i = higher; i >= lower; i--)
		{
			if (vout(i - 1) > vout(i))
			{
				swap(vout(i - 1), vout(i));
				swap(index1(i - 1), index1(i));
				swap(index2(i - 1), index2(i));
				lastswap = i;
			}
		}
		lower = lastswap + 1;

		// bubble up
		for (int i = lower; i <= higher; i++)
		{
			if (vout(i - 1) > vout(i))
			{
				swap(vout(i - 1), vout(i));
				swap(index1(i - 1), index1(i));
				swap(index2(i - 1), index2(i));
				lastswap = i;
			}
		}
		higher = lastswap - 1;

	} while (lower <= higher);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 * @return 
 */
double MathTools::median(Array<double, 1> v)
{
	MathTools::shakersort(v);

	int n = v.size();
	if (n % 2 == 0)
		return (v(n / 2) + v((n / 2) - 1)) / 2.;
	else
		return v((n - 1) / 2);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param v
 * @return 
 */
double MathTools::median(Array<float, 1> v)
{
	MathTools::shakersort(v);

	int n = v.size();
	if (n % 2 == 0)
		return (v(n / 2) + v((n / 2) - 1)) / 2.;
	else
		return v((n - 1) / 2);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param n
 * @param v
 * @return 
 */
double MathTools::closestValue(double n, Array<double, 1> v)
{
	double d0 = 1, d1 = 0;
	unsigned int i;
	for (i = 0; i < v.size() - 1 && d0 > d1; i++)
	{
		d0 = fabs(n - v(i));
		d1 = fabs(n - v(i + 1));
	}
	if (d0 <= d1)
		return v(i - 1);
	else
		return v(i);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param n
 * @param v
 * @param cv
 * @param index
 */
void MathTools::closestValue(double n, Array<double, 1> v, double &cv, int &index)
{
	double d0 = 1, d1 = 0;
	unsigned int i;
	for (i = 0; i < v.size() - 1 && d0 > d1; i++)
	{
		d0 = fabs(n - v(i));
		d1 = fabs(n - v(i + 1));
	}
	if (d0 <= d1)
	{
		cv = v(i - 1);
		index = i - 1;
	}
	else
	{
		cv = v(i);
		index = i;
	}
}

//==============================================================================




//==============================================================================
/**
 * 
 * @param x1
 * @param y1
 * @param x2
 * @param y2
 * @param x
 * @return 
 */
double MathTools::interpolate(double x1, double y1, double x2, double y2, double x)
{
	return (y2 - (x2 - x) * (y2 - y1) / (x2 - x1));

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xmask
 * @param x
 * @param y
 * @param temp
 */
void MathTools::interpolate(Array<double, 1> xmask, Array<double, 1> x, Array<double, 1> y, Array<double, 1> &temp)
{
	//#*****Alle Spektra werden linear interpoliert.
	//Maske der x und deren Anzahl sind xmask bzw. num_xmask
	//Zu interpolierende x und y-Werte und deren Anzahl sind x,y und num_data
	//interpolierte Werte und Anzahl sind x,y und num_ip

	temp.resize(xmask.size());
	double x1, x2, y1, y2;

	//sort everything
	MathTools::shakersort(xmask);
	MathTools::shakersort(x, y);
	//Search for values in the mask that are smaler than the first data value
	unsigned int ll = 0;
	while ((xmask(ll) <= x(0)) && (ll < xmask.size()))
	{
		ll++;
	};
	//set these values to the first function-value of the data
	for (unsigned int j = 0; j < ll; j++)
		temp(j) = y(0);

	unsigned int j, k = 0;
	int k0 = 0;

	//loop through the mask values
	for (j = ll; j < xmask.size(); j++)
	{
		//search for data values, between which a mask value lies
		for (k = k0; k < x.size() - 1; k++)
			if (x(k) <= xmask(j) && x(k + 1) > xmask(j))
				break;
		//if the mask extends beyond the data stop
		if (k == x.size() - 1)
			break;

		//save the momentary x-index to proceed at that point in the next loop to save time
		k0 = k;

		x1 = x(k);
		y1 = y(k);
		x2 = x(k + 1);
		y2 = y(k + 1);

		//interpolation
		temp(j) = y2 - (x2 - xmask(j)) * (y2 - y1) / (x2 - x1);
	}

	//copy nonexisting values at the end  of the data
	for (unsigned int i = j; i < xmask.size(); i++)
		temp(j) = y(y.size() - 1);
}

/*
//==============================================================================




//==============================================================================
 void MathTools::getRange ( QVector<double> &v1, double min, double max )
 {
 QVector<double> range;
 for ( int i = 0;i < v1.size();i++ )
 if ( v1[i] >= min && v1[i] <= max ) range.append ( v1[i] );

 v1 = range;
 }
//==============================================================================




//==============================================================================
 void MathTools::getRange ( QVector<double> &v1, QVector<double> &v2, double minv1, double maxv1 )
 {
 QVector<double> v1cut, v2cut;
 for ( int i = 0;i < v1.size();i++ )
 if ( v1[i] >= minv1 && v1[i] <= maxv1 )
 {
 v1cut.append ( v1[i] );
 v2cut.append ( v2[i] );
 }

 v1 = v1cut;
 v2 = v2cut;
 }
//==============================================================================




//==============================================================================
 void MathTools::compact ( QVector<double> f, QVector<double> a,  QVector<double> &fc, QVector<double> &ac )
 {
 fc.clear();
 ac.clear();

 if ( f.isEmpty() ) return;

 fc.append ( f[0] );
 ac.append ( a[0] );
 for ( int i = 1;i < f.size() - 1;i++ )
 if ( MathTools::signof ( a[i] - a[i-1] ) != MathTools::signof ( a[i+1] - a[i] ) )
 {
 fc.append ( f[i] );
 ac.append ( a[i] );
 }
 }
//==============================================================================




//==============================================================================
 bool MathTools::getLombScargle ( QVector<double> x, QVector<double> y, QVector<double> w, double nu_min, double nu_max, double dnu, QVector<double> &freq, QVector<double> &power, QVector<double> &amplitude, double &maxAmp, double &maxFreq, QProgressDialog &progress )
 {
 if ( fabs ( nu_min ) < 0.0001 ) nu_min = 0.0001;

 int numNu = int ( ( nu_max - nu_min ) / dnu );
 double nu_min2Pi = Constants::Pi2 * nu_min;
 double dnu2Pi = Constants::Pi2 * dnu;
 double numx_double = double ( x.size() );
 double sqrNumx = numx_double * numx_double;
 QVector<double> ss ( numNu, 0. ), ss2 ( numNu, 0. ), sc ( numNu, 0. ), sc2 ( numNu, 0. );
 double AF0, sinAF0, cosAF0, S20, C20, ADF, SDF, CDF, S2DF, C2DF, SSK, SS2K, SCK, SC2K, C0X, S0X, CTX, C2T;
 double pfacanzpkt = 4. / x.size();

 double wSum = 0.;
 for ( int j = 0;j < x.size();j++ )  wSum += w[j];
 double numxdWsum = x.size() / wSum;

 freq.resize ( numNu );
 power.resize ( numNu );
 amplitude.resize ( numNu );

 double mean = Statistics::getMean ( y, w );
 double yave;
 maxAmp = -999999.;

 progress.setRange ( 0, x.size() );
 progress.setModal ( true );
 progress.setMinimumDuration ( 0 );

 for ( int i = 0;i < x.size();i++ )
 {
 if ( progress.labelText() != "invisible" )
 {
 progress.setValue ( i );
 if ( progress.wasCanceled() ) return false;
 }

 yave = y[i] - mean;

 AF0 = fmod ( x[i] * nu_min2Pi, Constants::Pi2 );
 sinAF0 = sin ( AF0 );
 cosAF0 = cos ( AF0 );
 S20 = 2. * sinAF0 * cosAF0;
 C20 = cosAF0 * cosAF0 - sinAF0 * sinAF0;

 ADF = fmod ( x[i] * dnu2Pi, Constants::Pi2 );
 SDF = sin ( ADF );
 CDF = cos ( ADF );
 S2DF = 2. * SDF * CDF;
 C2DF = CDF * CDF - SDF * SDF;

 C0X = cosAF0 * w[i] * yave;
 S0X = sinAF0 * w[i] * yave;

 for ( int k = 0;k < numNu;k++ )
 {
 ss [k] += S0X;
 sc [k] += C0X;
 CTX = C0X;
 C0X = CTX * CDF - S0X * SDF;
 S0X = S0X * CDF + CTX * SDF;
 ss2 [k] += S20;
 sc2 [k] += C20;
 C2T = C20;
 C20 = C2T * C2DF - S20 * S2DF;
 S20 = S20 * C2DF + C2T * S2DF;
 }

 }

 double nu = nu_min;

 for ( int k = 0;k < numNu;k++ )
 {
 SSK  =  ss [k] * numxdWsum;
 SS2K = ss2 [k] * numxdWsum;
 SCK  =  sc [k] * numxdWsum;
 SC2K = sc2 [k] * numxdWsum;

 freq[k] = nu;
 if ( nu == 0 )
 {
 power[k] = 0.;
 amplitude[k] = 0.;
 }
 power[k] =  pfacanzpkt * ( SCK * SCK * ( numx_double - SC2K ) + SSK * SSK * ( numx_double + SC2K ) - 2. * SSK * SCK * SS2K ) /    ( sqrNumx - SC2K * SC2K - SS2K * SS2K );
 amplitude[k] =  sqrt ( power[k]  );

 if ( amplitude[k] > maxAmp )
 {
 maxAmp = amplitude[k];
 maxFreq = freq[k];
 }
 nu += dnu;
 }

 return true;
 }

 */

//==============================================================================




//==============================================================================
/**
 * 
 * @param xe
 * @param ze
 * @param xt
 * @param yt
 * @param zt
 * @param siwi
 */
void MathTools::getLOSAngle(double xe, double ze, double xt, double yt, double zt, double &siwi)
{
	double dummy;

	dummy = (xt * xe + zt * ze) / sqrt((xt * xt + yt * yt + zt * zt));
	siwi = acos(dummy);

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 * @param radius
 * @param theta
 * @param fi
 */
void MathTools::karth2Polar(double x, double y, double z, double &radius, double &theta, double &fi)
{
	//theta is declination! (not zenith angle)
	radius = sqrt(x * x + y * y + z * z);
	if (radius == 0.)
	{
		theta = fi = 0.;
		return;
	}

	if (y >= 0)
		fi = acos(x / (sqrt(x * x + y * y)));
	else
		fi = -acos(x / (sqrt(x * x + y * y)));
	theta = asin(z);
	if (fi < 0.)
		fi += Constants::Pi2;
	if (fi > Constants::Pi2)
		fi -= Constants::Pi2;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 * @param radius
 * @param theta
 * @param fi
 */
void MathTools::karth2Polar(float x, float y, float z, float &radius, float &theta, float &fi)
{
	//theta is declination! (not zenith angle)
	radius = sqrt(x * x + y * y + z * z);
	if (radius == 0.)
	{
		theta = fi = 0.;
		return;
	}

	if (y >= 0)
		fi = acos(x / (sqrt(x * x + y * y)));
	else
		fi = -acos(x / (sqrt(x * x + y * y)));
	theta = asin(z);
	if (fi < 0.)
		fi += Constants::Pi2;
	if (fi > Constants::Pi2)
		fi -= Constants::Pi2;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param radius
 * @param theta
 * @param fi
 * @param X
 * @param Y
 * @param Z
 */
void MathTools::polar2Karth(double radius, double theta, double fi, double &X, double &Y, double &Z)
{
	double rsint;
	rsint = radius * sin(theta);
	X = rsint * cos(fi);
	Y = rsint * sin(fi);
	Z = radius * cos(theta);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param radius
 * @param theta
 * @param fi
 * @param X
 * @param Y
 * @param Z
 */
void MathTools::polar2Karth(float radius, float theta, float fi, float &X, float &Y, float &Z)
{
	float rsint;
	rsint = radius * sin(theta);
	X = rsint * cos(fi);
	Y = rsint * sin(fi);
	Z = radius * cos(theta);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fitime
 * @param e
 * @param radius
 * @param theta
 * @param fi
 * @param X
 * @param Y
 * @param Z
 */
void MathTools::polarTrans2Karth(double fitime, double e, double radius, double theta, double fi, double &X, double &Y, double &Z)
{
	double rsint, costheta, sintheta, sinfi, cosfi;

	costheta = cos(theta) * cos(e) + sin(theta) * sin(e) * cos(fi);
	sintheta = sqrt(1 - costheta * costheta);
	sinfi = sin(fi) * sin(theta) / sintheta;
	cosfi = (cos(e) * sin(theta) * cos(fi) - sin(e) * cos(theta)) / sintheta;

	rsint = radius * sintheta;
	X = rsint * (cosfi * cos(fitime) - sinfi * sin(fitime));
	Y = rsint * (sinfi * cos(fitime) + cosfi * sin(fitime));
	Z = radius * costheta;
	//  if ((cosfi>0) && (sintheta>0)) {(*X)=-(*X);}

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xe
 * @param ze
 * @param X
 * @param Y
 * @param Z
 * @param A
 * @param B
 */
void MathTools::trans3D2D(double xe, double ze, double X, double Y, double Z, double &A, double &B)
{
	A = Y;
	B = (Z * xe - X * ze);
}

/*

 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 void MathTools::areaPolygon ( QVector<double> x, QVector<double> y, int AnzPts, double &Area )
 {
 int i;
 AnzPts--;
 Area = x[0] * ( y[1] - y[AnzPts] ) + x[AnzPts] * ( y[0] - y[AnzPts-1] );
 for ( i = 1;i < AnzPts;i++ )
 {
 Area += x[i] * ( y[i+1] - y[i-1] );
 }
 Area = 0.5 * fabs ( Area );
 }
 */
//==============================================================================




//==============================================================================
/**
 * 
 * @param x0
 * @param y0
 * @param z0
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param Xn
 * @param Yn
 * @param Zn
 */
void MathTools::crossProduct(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double &Xn,
		double &Yn, double &Zn)
{
	Xn = (y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0);
	Yn = (z1 - z0) * (x2 - x0) - (x1 - x0) * (z2 - z0);
	Zn = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x0
 * @param y0
 * @param z0
 * @param x1
 * @param y1
 * @param z1
 * @param Xn
 * @param Yn
 * @param Zn
 */
void MathTools::crossProduct(double x0, double y0, double z0, double x1, double y1, double z1, double &Xn, double &Yn, double &Zn)
{
	Xn = y0 * z1 - y1 * z0;
	Yn = -x0 * z1 + x1 * z0;
	Zn = x0 * y1 - x1 * y0;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x0
 * @param y0
 * @param z0
 * @param x1
 * @param y1
 * @param z1
 * @param Xn
 * @param Yn
 * @param Zn
 */
void MathTools::crossProduct(float x0, float y0, float z0, float x1, float y1, float z1, float &Xn, float &Yn, float &Zn)
{
	Xn = y0 * z1 - y1 * z0;
	Yn = -x0 * z1 + x1 * z0;
	Zn = x0 * y1 - x1 * y0;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 * @return 
 */
double MathTools::absoluteValue(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 * @return 
 */
double MathTools::absoluteValue(float x, float y, float z)
{
	return sqrt(x * x + y * y + z * z);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @return 
 */
double MathTools::absoluteValue(double x, double y)
{
	return sqrt(x * x + y * y);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 */
void MathTools::normalize(double &x, double &y, double &z)
{
	double abs = MathTools::absoluteValue(x, y, z);
	x /= abs;
	y /= abs;
	z /= abs;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x
 * @param y
 * @param z
 */
void MathTools::normalize(float &x, float &y, float &z)
{
	double abs = MathTools::absoluteValue(x, y, z);
	x /= abs;
	y /= abs;
	z /= abs;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @param x3
 * @param y3
 * @param z3
 */
void MathTools::vectorDiff(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3)
{
	x3 = x1 - x2;
	y3 = y1 - y2;
	z3 = z1 - z2;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param funcvalue
 * @param funcvalue_dx
 * @param dx
 * @return 
 */
double MathTools::derivative(double funcvalue, double funcvalue_dx, double dx)
{
	return (funcvalue_dx - funcvalue) / dx;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param theta
 * @param fi
 * @param thetaout
 * @param fiout
 */
void MathTools::correctSphAnglesOutOfLimit(double theta, double fi, double &thetaout, double &fiout)
{
	if (theta < 0)
	{
		thetaout = -theta;
		fiout = fi + Constants::Pi;
	}
	else if (theta > Constants::Pi)
	{
		thetaout = Constants::Pi2 - theta;
		fiout = fi + Constants::Pi;
	}
}
//==============================================================================




//==============================================================================
/**
 * This function computes the cross product of 2 vectors. The 2 vectors are defined in spherical coordinates. 
 * The origin of the 2 vector is defined by vector1(r1,t1,f1). The ends of the 2 vectors are defined by vector2 and 3. 
 * We convert the sph. coord. into carth. coord. The output are carthesian coordinates.
 * @param r1
 * @param t1
 * @param f1
 * @param r2
 * @param t2
 * @param f2
 * @param r3
 * @param t3
 * @param f3
 * @param x
 * @param y
 * @param z
 */
void MathTools::crossProductSph(double r1, double t1, double f1, double r2, double t2, double f2, double r3, double t3, double f3,
		double &x, double &y, double &z)
//
{
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;

	polar2Karth(r1, t1, f1, x1, y1, z1);
	polar2Karth(r2, t2, f2, x2, y2, z2);
	polar2Karth(r3, t3, f3, x3, y3, z3);

	crossProduct(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param x1
 * @param y1
 * @param x2
 * @param y2
 * @return 
 */
double MathTools::distance(double x1, double y1, double x2, double y2)
{
	return (sqrt(pow(x2 - x1, 2.) + pow(y2 - y1, 2.)));

}
//==============================================================================




//==============================================================================
/**
 * compute the area of a circle segment. r=radius of circle, phi=segment angle (in radial direction) (Netz, Formeln der Mathematik)
 * @param r
 * @param phi
 * @return 
 */
double MathTools::getCircleSegmentArea(double r, double phi)
{
	//double phi = 2. * acos ( 1 - h / r );  //segment angle
	return (r * r / 2. * (phi - sin(phi))); //= area
}
//==============================================================================




//==============================================================================
/**
 * compute the visible are of a disk with radius r1 that is eclipsed by another disk with radius r2
 * IMPORTANT: r2 is the smaller circle!
 * @param x1
 * @param y1
 * @param r1
 * @param x2
 * @param y2
 * @param r2
 * @return 
 */
double MathTools::getCircleOverlapArea(double x1, double y1, double r1, double x2, double y2, double r2)
{
	double darea = 0; //overlap area of the disks
	double dist = MathTools::distance(x1, y1, x2, y2);
	double phi1, phi2;
	if (dist < (r1 + r2) && dist > (r1 - r2)) //condition: partial overlap of disks
	{
		phi1 = 2. * MathTools::getAngleTriangle(r2, r1, dist);
		phi2 = 2. * MathTools::getAngleTriangle(r1, r2, dist);
		darea = MathTools::getCircleSegmentArea(r1, phi1) + MathTools::getCircleSegmentArea(r2, phi2);
	}
	else if (dist <= (r1 - r2))
		darea = r2 * r2 * Constants::Pi;
	return darea;
}

//==============================================================================




//==============================================================================
/**
 * compute the angle in a triangle with sides a,b,c that is opposite to side a
 * @param a
 * @param b
 * @param c
 * @return 
 */
double MathTools::getAngleTriangle(double a, double b, double c)
{
	return (acos((-a * a + b * b + c * c) / (2. * b * c)));
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param fi1
 * @param theta1
 * @param fi2
 * @param theta2
 * @return 
 */
double MathTools::getAngularDistanceSphere(double fi1, double theta1, double fi2, double theta2)
{
	return acos(sin(theta1) * sin(theta2) + cos(theta1) * cos(theta2) * cos(fi2 - fi1));

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param X0
 * @param Y0
 * @param eTh
 * @param eFi
 * @param sTh
 * @param sFi
 */
void MathTools::trans2DPolar(double X0, double Y0, double eTh, double eFi, double &sTh, double &sFi)
{
	//compute the polar coordinates
	double X, Y, Z, r;

	//when eTh==0 then the calculation is not defined (since X and Y are derived by divisions through sin(eTh)
	if (eTh == 0)
		eTh = 1E-10;

	double sineTh = sin(eTh);
	double coseTh = cos(eTh);
	double sineFi = sin(eFi);
	double coseFi = cos(eFi);

	double xP1 = sineTh * X0 * sineFi - coseFi * Y0 + coseFi * pow2(coseTh) * Y0;
	double xP2 = coseFi * coseTh * sqrt(-pow2(sineTh) * (Y0 * Y0 + X0 * X0 - 1));
	double yP1 = -sineTh * coseFi * X0 + pow2(coseTh) * sineFi * Y0 - Y0 * sineFi;
	double yP2 = coseTh * sineFi * sqrt(-pow2(sineTh) * (Y0 * Y0 + X0 * X0 - 1));
	double zP1 = coseTh * Y0;
	double zP2 = sqrt(-pow2(sineTh) * (Y0 * Y0 + X0 * X0 - 1));

	X = (xP1 + xP2) / sineTh;
	Y = (yP1 + yP2) / sineTh;
	Z = zP1 + zP2;

	karth2Polar(X, Y, Z, r, sTh, sFi);

	if (getAngularDistanceSphere(eFi, eTh, sFi, sTh) > Constants::Pid2)
	{
		X = (xP1 - xP2) / sineTh;
		Y = (yP1 - yP2) / sineTh;
		Z = zP1 - zP2;

		karth2Polar(X, Y, Z, r, sTh, sFi);
	}

}
//==============================================================================




//==============================================================================
/**
 * returns the normal distance of the point from the line
 * from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
 * @param xp point
 * @param yp point
 * @param zp point
 * @param ux line vector
 * @param uy line vector
 * @param uz line vector
 * @return distance
 */
double MathTools::distancePointLine(double xp, double yp, double zp, double ux, double uy, double uz)
{

	double xn, yn, zn;
	crossProduct(xp, yp, zp, xp - ux, yp - uy, zp - uz, xn, yn, zn);
	double d = absoluteValue(xn, yn, zn) / absoluteValue(ux, uy, uz);

	return d;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xp
 * @param yp
 * @param zp
 * @param ux
 * @param uy
 * @param uz
 * @param theta
 */
void MathTools::rotation(double &xp, double &yp, double &zp, double ux, double uy, double uz, double theta)
{
	// xp, yp, zp;  //vector rotated by angle theta around axis (ux,uy,uz)

	double c = cos(theta);
	double s = sin(theta);

	normalize(ux, uy, uz);

	//from http://en.wikipedia.org/wiki/Rotation_matrix
	double xd = xp * (ux * ux + (1 - ux * ux) * c) + yp * (ux * uy * (1 - c) - uz * s) + zp * (ux * uz * (1 - c) + uy * s);
	double yd = xp * (ux * uy * (1 - c) + uz * s) + yp * (uy * uy + (1 - uy * uy) * c) + zp * (uy * uz * (1 - c) - ux * s);
	double zd = xp * (ux * uz * (1 - c) - uy * s) + yp * (uy * uz * (1 - c) + ux * s) + zp * (uz * uz + (1 - uz * uz) * c);

	xp = xd;
	yp = yd;
	zp = zd;
}
//==============================================================================




//==============================================================================
/**
 * vector rotated by angle theta around axis (ux,uy,uz)
 * @param xp
 * @param yp
 * @param zp
 * @param ux
 * @param uy
 * @param uz
 * @param theta rotation angle
 */
void MathTools::rotation(float &xp, float &yp, float &zp, float ux, float uy, float uz, float theta)
{

	float c = cos(theta);
	float s = sin(theta);

	normalize(ux, uy, uz);

	//from http://en.wikipedia.org/wiki/Rotation_matrix
	float xd = xp * (ux * ux + (1 - ux * ux) * c) + yp * (ux * uy * (1 - c) - uz * s) + zp * (ux * uz * (1 - c) + uy * s);
	float yd = xp * (ux * uy * (1 - c) + uz * s) + yp * (uy * uy + (1 - uy * uy) * c) + zp * (uy * uz * (1 - c) - ux * s);
	float zd = xp * (ux * uz * (1 - c) - uy * s) + yp * (uy * uz * (1 - c) + ux * s) + zp * (uz * uz + (1 - uz * uz) * c);

	xp = xd;
	yp = yd;
	zp = zd;
}

//==============================================================================




//==============================================================================
/**
 * 
 * @param xp point to be projected
 * @param yp point to be projected
 * @param zp point to be projected
 * @param a, b, c normal vector of plane (ax+by+cz=0) 
 * @param xd normal projected point in plane
 * @param yd normal projected point in plane
 * @param zd normal projected point in plane
 */
void MathTools::normalProjectionPointOnPlane(double xp, double yp, double zp, double a, double b, double c, double &xd, double &yd,
		double &zd)
{
	double invabcsq = 1. / (a * a + b * b + c * c);
	double ab = a * b;
	double ac = a * c;
	double bc = b * c;

	//from http://www.vsmp.ch/de/bulletins/no98/keller.pdf
	xd = invabcsq * (xp * (b * b + c * c) - yp * ab - zp * ac);
	yd = invabcsq * (-xp * ab + yp * (c * c + a * a) - zp * bc);
	zd = invabcsq * (-xp * ac - yp * bc + zp * (a * a + b * b));

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xp point to be projected
 * @param yp point to be projected
 * @param zp point to be projected
 * @param a, b, c normal vector of plane (ax+by+cz=0) 
 * @param xd normal projected point in plane
 * @param yd normal projected point in plane
 * @param zd normal projected point in plane
 */
void MathTools::normalProjectionPointOnPlane(float xp, float yp, float zp, float a, float b, float c, float &xd, float &yd, float &zd)
{
	double invabcsq = 1. / (a * a + b * b + c * c);
	double ab = a * b;
	double ac = a * c;
	double bc = b * c;

	//from http://www.vsmp.ch/de/bulletins/no98/keller.pdf
	xd = invabcsq * (xp * (b * b + c * c) - yp * ab - zp * ac);
	yd = invabcsq * (-xp * ab + yp * (c * c + a * a) - zp * bc);
	zd = invabcsq * (-xp * ac - yp * bc + zp * (a * a + b * b));

}
//==============================================================================






//==============================================================================
/**
 * 
 * @param xp
 * @param yp
 * @param zp
 * @param ix
 * @param iy
 * @param iz
 * @param jx
 * @param jy
 * @param jz
 * @param kx
 * @param ky
 * @param kz
 * @param xt
 * @param yt
 * @param zt
 */
void MathTools::coordinateTransformation(double xp, double yp, double zp, double ix, double iy, double iz, double jx, double jy, double jz,
		double kx, double ky, double kz, double &xt, double &yt, double &zt)
{
	// xp, yp, zp;  point to be transformed
	// i, j, k;   normalized coordinates of axes of coordinate system of point
	// xt, yt, zt;  transformed point

	//from http://www.kwon3d.com/theory/transform/transform.html
	xt = xp * ix + yp * iy + zp * iz;
	yt = xp * jx + yp * jy + zp * jz;
	zt = xp * kx + yp * ky + zp * kz;

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param ra Right Ascension of object in radians
 * @param dec Declination in radians
 * @param ra0 R.A. of (x,y)=(0,0) in radians
 * @param dec0 Dec of (x,y)=(0,0) in radians
 * @param gamma Rotation angle of field in rad
 * @param x x-coordinate (on CCD or plane)
 * @param y y-coordinate (on CCD or plane)
 */
void MathTools::getGnomonicProjection(double ra, double dec, double ra0, double dec0, double gamma, double &x, double &y)
{

	//dx...    x Offset in projection plane
	//dy...    y Offset in projection plane


	double dummy = cos(dec0) * cos(dec) * cos(ra - ra0) + sin(dec0) * sin(dec);
	double X = (cos(dec) * sin(ra - ra0)) / dummy;
	double Y = -(sin(dec0) * cos(dec) * cos(ra - ra0) - cos(dec0) * sin(dec)) / dummy;

	x = -X * cos(gamma) - Y * sin(gamma);// - dx;
	y = -X * sin(gamma) + Y * cos(gamma);// - dy;
    
    
    }
//==============================================================================




//==============================================================================
/**
 * 
 * @param ra Right Ascension of object in radians
 * @param dec Declination in radians
 * @param ra0 R.A. of (x,y)=(0,0) in radians
 * @param dec0 Dec of (x,y)=(0,0) in radians
 * @param gamma Rotation angle of field in rad
 * @param x x-coordinate (on CCD or plane)
 * @param y y-coordinate (on CCD or plane)
 */
void MathTools::getGnomonicProjection(float ra, float dec, float ra0, float dec0, float gamma, float &x, float &y)
{
	//dx...    x Offset in projection plane
	//dy...    y Offset in projection plane


	float dummy = cos(dec0) * cos(dec) * cos(ra - ra0) + sin(dec0) * sin(dec);
	float X = (cos(dec) * sin(ra - ra0)) / dummy;
	float Y = -(sin(dec0) * cos(dec) * cos(ra - ra0) - cos(dec0) * sin(dec)) / dummy;

	x = -X * cos(gamma) - Y * sin(gamma);// - dx;
	y = -X * sin(gamma) + Y * cos(gamma);// - dy;


}
//==============================================================================




//==============================================================================
/**
 * 
 * @param ra Right Ascension of object in radians
 * @param dec Declination in radians
 * @param ra0 R.A. of (x,y)=(0,0) in radians
 * @param dec0 Dec of (x,y)=(0,0) in radians
 * @param gamma Rotation angle of field in rad
 * @param x x-coordinate (on CCD or plane)
 * @param y y-coordinate (on CCD or plane)
 */
void MathTools::getInverseGnomonicProjection(double x, double y, double ra0, double dec0, double gamma, double &ra, double &dec)
{
	//x += dx;
	//y += dy;
    	//dx...    x Offset in projection plane
	//dy...    y Offset in projection plane
	double X = x * cos(gamma) - y * sin(gamma);
	double Y = x * sin(gamma) + y * cos(gamma);

	ra = ra0 + atan2(-X, cos(dec0) - Y * sin(dec0));
	dec = asin((sin(dec0) + Y * cos(dec0)) / sqrt(1 + X * X + Y * Y));

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param ra Right Ascension of object in radians
 * @param dec Declination in radians
 * @param ra0 R.A. of (x,y)=(0,0) in radians
 * @param dec0 Dec of (x,y)=(0,0) in radians
 * @param gamma Rotation angle of field in rad
 * @param x x-coordinate (on CCD or plane)
 * @param y y-coordinate (on CCD or plane)
 */
void MathTools::getInverseGnomonicProjection(float x, float y, float ra0, float dec0, float gamma, float &ra, float &dec)
{
	//dx...    x Offset in projection plane
	//dy...    y Offset in projection plane
	float X = x * cos(gamma) - y * sin(gamma);
	float Y = x * sin(gamma) + y * cos(gamma);

	ra = ra0 + atan2(-X, cos(dec0) - Y * sin(dec0));
	dec = asin((sin(dec0) + Y * cos(dec0)) / sqrt(1 + X * X + Y * Y));

}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xIn
 * @param yIn
 * @param x0
 * @param y0
 * @param angle
 * @param xOut
 * @param yOut
 */
void MathTools::rotate(double xIn, double yIn, double x0, double y0, double angle, double &xOut, double &yOut)
{
	xOut = cos(angle) * (xIn - x0) - sin(angle) * (yIn - y0) + x0;
	yOut = sin(angle) * (xIn - x0) + cos(angle) * (yIn - y0) + y0;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param xIn
 * @param yIn
 * @param x0
 * @param y0
 * @param angle
 * @param xOut
 * @param yOut
 */
void MathTools::rotate(float xIn, float yIn, float x0, float y0, float angle, float &xOut, float &yOut)
{
	//rotate the coordinates (xIn, yIn) around the point (x0,y0) by the angle angle
	xOut = cos(angle) * (xIn - x0) - sin(angle) * (yIn - y0) + x0;
	yOut = sin(angle) * (xIn - x0) + cos(angle) * (yIn - y0) + y0;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param posIn x,y-coordinates of all stars on the CCD
 * @param x0 Offset of the rotation axis of the jitter roll relative to (0,0) on CCD
 * @param y0 Offset of the rotation axis of the jitter roll relative to (0,0) on CCD
 * @param yaw Offset of yaw in decimal pixels
 * @param pitch Offset of pitch in decimal pixels
 * @param roll Roll in radians
 * @param posOut Shifted x,y-coordinates of all stars of CCD
 */
void MathTools::rotate(Array<double, 2> posIn, double x0, double y0, double yaw, double pitch, double roll, Array<double, 2> &posOut)
{
	posOut(Range::all(), 0) = x0 + yaw + cos(roll) * (posIn(Range::all(), 0) - x0) - sin(roll) * (posIn(Range::all(), 1) - y0);
	posOut(Range::all(), 1) = y0 + pitch + sin(roll) * (posIn(Range::all(), 0) - x0) + cos(roll) * (posIn(Range::all(), 1) - y0);
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param posIn
 * @param x0
 * @param y0
 * @param yaw
 * @param roll
 * @param posOut
 */
void MathTools::rotate(Array<double, 2> posIn, double x0, double y0, double yaw, double roll, Array<double, 1> &posOut)
{
	posOut = yaw + cos(roll) * (posIn(Range::all(), 0) - x0) - sin(roll) * (posIn(Range::all(), 1) - y0) + x0;
}
//==============================================================================




//==============================================================================
/**
 * 
 * @param mapIn
 * @param x0
 * @param y0
 * @param angle
 * @param mapOut
 */
void MathTools::rotateMap(Array<double, 2> mapIn, double x0, double y0, double angle, Array<double, 2> &mapOut)
{
	//mapOut = 0.;
	double step = 0.5;
	//qDebug ( "void Jitter::rotateMap ( Array<int, 2> mapIn, Array<int, 2> mapOut )" );
	//qDebug() << "start rotation" << std::endl;
	double xrot, yrot;
	int xnew, ynew;
//	int counter = 0;
	for (double x = 0; x < mapIn.rows(); x++)
		for (double y = 0; y < mapIn.cols(); y++)
			for (double dx = 0; dx < 1; dx += step)
				for (double dy = 0; dy < 1; dy += step)
				{
					MathTools::rotate(x + dx, y + dy, x0, y0, angle, xrot, yrot);
					xnew = int(floor(xrot));
					ynew = int(floor(yrot));
					if (xnew >= 0 && xnew < mapIn.rows() && ynew >= 0 && ynew < mapIn.cols())
						mapOut(xnew, ynew) += mapIn(int(x), int(y));
				}

}
//==============================================================================




//==============================================================================
/**
 * Translation of a point (x,y) by (dx,dy).
 * Linear coordinate transformation of (x, y) with the shift of (dx, dy)
 * @param x
 * @param y
 * @param dx
 * @param dy
 */
void MathTools::translate(double &x, double &y, double dx, double dy)
{
	
	x -= dx;
	y -= dy;
}
//==============================================================================




//==============================================================================
/**
 * Translation of a point (x,y) by (dx,dy).
 * Linear coordinate transformation of (x, y) with the shift of (dx, dy)
 * @param x
 * @param y
 * @param dx
 * @param dy
 */
void MathTools::translate(float &x, float &y, float dx, float dy)
{
	
	x -= dx;
	y -= dy;
}
//==============================================================================




//==============================================================================
/**
 * computes the angle between two vectors in 3D
 * @param ux
 * @param uy
 * @param uz
 * @param vx
 * @param vy
 * @param vz
 * @return 
 */
double MathTools::getAngle(double ux, double uy, double uz, double vx, double vy, double vz)
{
	
	double dummy = MathTools::getInnerProduct(ux, uy, uz, vx, vy, vz) / (MathTools::absoluteValue(ux, uy, uz) * MathTools::absoluteValue(
			vx, vy, vz));

	if (dummy >= 1)
		return 0;
	else
		return acos(dummy);

}
//==============================================================================




//==============================================================================
/**
 * compute the angle in the plane between the two vectors (ux, uy), (vx, vy) in 2D
 * @param ux
 * @param uy
 * @param vx
 * @param vy
 * @return 
 */
double MathTools::getAngle(double ux, double uy, double vx, double vy)
{
	
	return acos(MathTools::getInnerProduct(ux, uy, vx, vy) / (MathTools::absoluteValue(ux, uy) * MathTools::absoluteValue(vx, vy)));
}
//==============================================================================




//==============================================================================
/**
 * Compute the inner product of two 3D-vectors.
 * @param ux
 * @param uy
 * @param uz
 * @param vx
 * @param vy
 * @param vz
 * @return 
 */
double MathTools::getInnerProduct(double ux, double uy, double uz, double vx, double vy, double vz)
{
	return ux * vx + uy * vy + uz * vz;
}
//==============================================================================




//==============================================================================
/**
 * Compute the inner product of two 2D-vectors.
 * @param ux
 * @param uy
 * @param vx
 * @param vy
 * @return 
 */
double MathTools::getInnerProduct(double ux, double uy, double vx, double vy)
{
	return ux * vx + uy * vy;
}
//==============================================================================




