///////////////////////////////////////////////////////////
// Tabfunction.h: Declaration AND implementation of the Tabfunction class
//                which allows to handle tabulated functions via automatic
//                interpolation between the tabulated values. Due to the
//                way many compilers handle templates, I need to include
//                also all the definitions...
//
//  Created on:      23-Oct-2012 1:59:58 PM
//  Original author: Joris De Ridder
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



#ifndef TABFUNCTION_H    // Include this file only once
#define TABFUNCTION_H




enum Extrap_Method {Linear_Extrapolation};
enum Interp_Method {Linear_Interpolation, Spline_Interpolation};



template <class Type>
class Tabfunction
{
 public:

  Tabfunction();
  ~Tabfunction();

  void Init(Type &a, Type &b, unsigned long N);
  void Set_Extrapolation_Method(Extrap_Method method);
  void Set_Interpolation_Method(Interp_Method method);

  double Begin();
  double End();

  double operator()(const double xvalue);

  double Integrate(double lower, double upper);
  void Set_Accuracy(const double accur);

 private:

  Extrap_Method extrap_method;
  Interp_Method interp_method;

  double accuracy;                // fraction accuracy for convergence

  unsigned long Nvalues;          // Number of tabulated values
  double *x;                      // The abscissa values
  double *y;                      // The ordinate values




  // Static variables used for spline interpolation

  double *DDy;                    // Second derivatives 
  double *temp;                   // temporary storage


  void init_spline();
  double spline_evaluate(const double xvalue);
  double rational_evaluate(const double xvalue);
  double linear_evaluate(const double xvalue);


};












//--------------------------------
  template <class Type>
  Tabfunction<Type>::Tabfunction()
//--------------------------------

// PURPOSE: default constructor...
// 
// REMARKS: . I initialise the array pointers to 0, so that I can harmless
//            delete [] them.
//
    : extrap_method(Linear_Extrapolation), interp_method(Spline_Interpolation),
    accuracy(1.0e-5), Nvalues(0), x(0), y(0), DDy(0), temp(0)
{

} // end Tabfunction()











//---------------------------------
  template <class Type>
  Tabfunction<Type>::~Tabfunction()
//---------------------------------

// PURPOSE: default destructor
{
  delete [] x;
  delete [] y;
  delete [] DDy;
  delete [] temp;


} // end ~Tabfunction()











//---------------------------------------------------------------------------
  template <class Type>
  void Tabfunction<Type>::Init(Type &a, Type &b, unsigned long N)
//---------------------------------------------------------------------------

// PURPOSE: Copy the user given arrays into internal arrays, and
//          initialize the interpolation procedure.
//
// INPUT: a: array with abscissa values of the tabulated function
//        b: corresponding array with ordinate values of the tabulated function
//        N: only elements 0..N-1 of a[] and b[] will be copied
//           to the internal arrays.
// 
// REMARK: . The abscissa array a must be sorted (ascending), and must be 
//           free from duplicate values.
//         . I use a template so that arrays, vectors and valarrays can
//           be used.
// 
{
  // Check whether the given size of the arrays is larger then zero
  // If so, allocate the necessary memory, if not complain.

  if (N > 0)
    {
      // First de-allocate memory. Necessary, if the user wants to
      // use Init multiple times

      Nvalues = N;
      delete [] x; 
      delete [] y;
      delete [] DDy;
      delete [] temp;

      // Then allocate the necessary memory.

      x = new double[N];
      y = new double[N];
      DDy = new double[N];
      temp = new double[N];
    }
  else
    {
      cerr << "\nError (Tabfunction::Init()): size <= 0" << endl;
      exit(1);
    }

  // Copy the user given arrays into internal arrays

  for (unsigned long i = 0; i < N; i++)
    {
      x[i] = double(a[i]);
      y[i] = double(b[i]);
    }


  // Check whether the abscissa array is sorted in an ascending way

  for (unsigned long i = 1; i < N; i++)
    {
      if (x[i] <= x[i-1])
	{
	  cerr << "\nError (Tabfunction::Init()):  x-values not sorted\n";
	  exit(1);
	}
    }

  // Initialize the spline

  init_spline();

} // end Init()














//----------------------------------------------------------------------
  template <class Type>
  void Tabfunction<Type>::Set_Extrapolation_Method(Extrap_Method method)
//----------------------------------------------------------------------

// PURPOSE: Decide how you wish to compute values out of the x-range,
//          by extrapolation.
{
  if (method == Linear_Extrapolation) 
    {
      extrap_method = method;
    }
  else
    {
      cerr << "\nError (Tabfunction::Set_Extrapolation_Method()) : "
           << "Illegal argument." << endl;
      exit(1);
    }
}















//----------------------------------------------------------------------
  template <class Type>
  void Tabfunction<Type>::Set_Interpolation_Method(Interp_Method method)
//----------------------------------------------------------------------

// PURPOSE: Decide how you wish to compute values inside the x-range,
//          by interpolation.
{
  if (   (method == Linear_Interpolation) 
      || (method == Spline_Interpolation))
    {
      interp_method = method;
    }
  else
    {
      cerr << "\nError (Tabfunction::Set_Interpolation_Method()):  "
           << "Illegal argument." << endl;
      exit(1);
    }
}

















//---------------------------------------------------------------
  template <class Type>
  double Tabfunction<Type>::operator()(const double xvalue)
//---------------------------------------------------------------

// PURPOSE: Evaluate the function by automatic interpolation in the
//          user given tabulated values.
//
{
  //  First check whether the class has already been initialized.

  if (Nvalues == 0)
  {
      cerr << "\nError ( Tabfunction::()): First give tabulated values.\n";
      exit(1);
  }

  // Check whether the given value is in its appropriate boundaries
  // and inter/extra-polate accordingly

  if ((xvalue >= x[0]) && (xvalue <= x[Nvalues-1]))
  {
      // Interpolate 

      if (interp_method == Linear_Interpolation)
      {
          return(linear_evaluate(xvalue));
      }

      if (interp_method == Spline_Interpolation)
      {
          return(spline_evaluate(xvalue));
      }
  }
  else
  {
      // Extrapolate

      if (extrap_method == Linear_Extrapolation)
      {
          return(linear_evaluate(xvalue));
      }
      
  
  }
  
  exit(0);

} // end operator () overloading
















//-------------------------------------
  template <class Type>
  void Tabfunction<Type>::init_spline()
//-------------------------------------

// PURPOSE: Initialize the necessary values for natural spline interpolation.
//          That is, compute the second derivatives of the interpolating 
//          function at the tabulated points. By "natural" spline we mean
//          that the second derivatives at the boundary values are assumed
//          to be zero.
//
{
  double p, qn, sig, un;

  // Set the lower boundary to be "natural", i.e. zero second derivative

  DDy[0] = temp[0] = 0.0;

  // Now the decomposition loop of the tridiagonal algorithm. DDy and u
  // are used for temporary storage of the decomposed factors.

  for (unsigned long i=1; i <= Nvalues-2; i++)
    { 
      sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
      p = sig * DDy[i-1] + 2.0;
      DDy[i] = (sig-1.0) / p;
      temp[i] =   (y[i+1]-y[i]) / (x[i+1]-x[i]) 
                - (y[i]-y[i-1]) / (x[i]-x[i-1]);
      temp[i] = (6.0 * temp[i] / (x[i+1]-x[i-1]) - sig * temp[i-1]) / p;
    }

  // Also set the upper boundary to be the "natural" boundary

  qn = un = 0.0; 
  DDy[Nvalues-1]=(un - qn * temp[Nvalues-2]) / (qn * DDy[Nvalues-2] + 1.0);

  // Now the backsubstitution loop of the tridiagonal algorithm
  // Pitfall: don't declare k unsigned, otherwise it will never reach
  // the value -1 to stop the loop.

  for (long k = Nvalues-2; k >= 0; k--)
    { 
      DDy[k] = DDy[k] * DDy[k+1] + temp[k]; 
    }


} // end init_spline()













//--------------------------------------------------------------------
  template <class Type>
  double Tabfunction<Type>::spline_evaluate(const double xvalue)
//--------------------------------------------------------------------

// PURPOSE: computes the natural cubic spline interpolation, given
//          the tabulated values of the function, and given the array
//          of second derivatives computed by init_spline().
//          "Natural" here means that the second derivatives of the function
//          at the tabulated boundary values are assumed to be zero. 
//
// SOURCE: Numerical Recipes
//
{
  unsigned long klo, khi, k;
  double h,b,a,yvalue;

  // First find the right place in the table by means of bisection.

  klo = 0; khi = Nvalues-1;
  while (khi-klo > 1)
    { 
      k = (khi+klo) >> 1;
      if (x[k] > xvalue) 
	{ 
	  khi = k;
	}
      else 
	{
	  klo = k;
	}
    }

  // Define a few handy abbriviations

  h = x[khi]-x[klo];
  a = (x[khi]-xvalue)/h;
  b = (xvalue-x[klo])/h;

  // Evaluate cubic spline

  yvalue =    a*y[klo] + b*y[khi] 
           + ((a*a*a-a)*DDy[klo]+(b*b*b-b)*DDy[khi]) * (h*h)/6.0;

  return(yvalue);

} // end spline_evaluate()














//--------------------------------------------------------------------
  template <class Type>
  double Tabfunction<Type>::linear_evaluate(const double xvalue)
//--------------------------------------------------------------------

// PURPOSE: perform linear inter- or extrapolation of a tabulated
//          function
{
  long i,j;      // Indices of the two points, defining the linear relation
  long m;        // Index of middle point

  // Out of range, or on the upper border

  if (xvalue >= x[Nvalues-1])
    {
      i = Nvalues-2;
      j = Nvalues-1;
      return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));
    }

  // Out of range, or on the lower border

  if (xvalue <= x[0])
    {
      i = 0;
      j = 1;
      return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));
    }

  // In tabulated range. First locate neighbours by bisection.

  i = 0;                 // We already checked the lower and upper
  j = Nvalues-1;         // borders.

  while (j - i > 1)
    {
      m = (j + i) >> 1;          // middle point

      if (xvalue >= x[m])
	{
	  i = m;
	}
      else
	{
	  j = m;
	}
    }

  return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));

} // end linear_evaluate()













//---------------------------------------------------------------------
  template <class Type>
  double Tabfunction<Type>::Integrate(double lower, double upper)
//---------------------------------------------------------------------

// PURPOSE: Integrate the tabulated function. The integrand is computed
//          via inter/extrapolation.
//
// INPUT: lower: lower integration boundary
//        upper: upper integration boundary
//
// OUTPUT: the integral from lower to upper.
//
// REMARKS: . source: Numerical Recipes, pg. 137.
//
{
  const int JMAX = 20;
  double x,tnm,sum,del,olds;
  double s = 0.0;
  int it,j,k;

  // Init old s to a number that is unlikely to the average of the function
  // at its endpoints

  olds = -1.0e30;

  // Iterate with finer trapezium coverages until converged

  for (j = 1; j <= JMAX; j++)
    {
      if (j == 1)
	{
	  s = 0.5 * (upper - lower) * ((*this)(lower) + (*this)(upper));
	}
      else
	{
	  for (it = 1, k = 1; k < j-1; k++) it <<= 1;
	  tnm = it;
	  del = (upper - lower)/tnm;
	  x = lower + 0.5 * del;
	  for (sum = 0.0, k=1; k <= it; k++, x+= del) sum += (*this)(x);
	  s = 0.5 * (s + (upper-lower)*sum / tnm);
	}

      // Avoid spurious early convergence, that's why the j > 5

      if (j > 5)
	{
	  if (fabs(s-olds) < accuracy * fabs(olds) ||
              (s == 0.0 && olds == 0.0)) return s;
	}
      
      olds = s;
    }


  cerr << "\nError (Tabfunction::Integrate):  no convergence." << endl;
  exit(1);

  return(0.0);   // never get here.

} // end Integrate()













//--------------------------------------------------------
  template <class Type>
  void Tabfunction<Type>::Set_Accuracy(const double accur)
//--------------------------------------------------------

// PURPOSE: set the fractional accuracy which determines when a convergence
//          procedure (e.g. to integrate) has to stop.
//
{
  if (accur > 0.0)
    {
      accuracy = accur;
    }
  else
    {
      cerr << "\nError (Tabfunction::Set_Accuracy()):  accuracy <= 0" << endl;
      exit(1);
    }

} // end Set_Accuracy()















//---------------------------------
  template <class Type>
  double Tabfunction<Type>::Begin()
//---------------------------------

// PURPOSE: return the abscissa value of the point with the lowest
//          abscissa value.
//
{
  //  First check whether the class has already been initialized.

  if (Nvalues == 0)
    {
      cerr << "\nError (Tabfunction::Begin()):  First give tabulated values.\n";
      exit(1);
    } 
  else
    {
      return(x[0]);
    }

} // end Begin()















//-------------------------------
  template <class Type>
  double Tabfunction<Type>::End()
//-------------------------------

// PURPOSE: return the abscissa value of the tabulated point with the greatest
//          abscissa value.
//
{
  //  First check whether the class has already been initialized.

  if (Nvalues == 0)
    {
      cerr << "\nError (Tabfunction::End()):  First give tabulated values.\n";
      exit(1);
    }
  else
    {
      return(x[Nvalues-1]);
    }

} // end Begin()






#endif









