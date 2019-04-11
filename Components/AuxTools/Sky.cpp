///////////////////////////////////////////////////////////
// Sky.cpp: C++ code to define a class that is able to return many properties
//         of the sky and objects in the sky like the sun.
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





#include "Sky.h"
#include "Skydata.h"


//==============================================================================
/**
 * Constructor method
 */
Sky::Sky()
      : obliquity ( 0.409087723 )
{

} // end Sky()
//==============================================================================
//==============================================================================
/**
 * Destructor method
 */
Sky::~Sky(){} // end ~Sky()
//==============================================================================






//==============================================================================
//Functions:
//==============================================================================
/**
 * Given a wavelength, this function returns the solar radiant flux
 * at this wavelength, measured above the atmosphere of the earth.
 * REMARKS: . This function uses the Wehrli (1985) solar irradiance table
 *            plus linear interpolation. Note that the units of the tabulated
 *            data radiant fluxes are J s^{-1} m^{-2} (nm)^{-1}, where nm
 *            is nanometer (unit of wavelength) while the units of the radiant
 *            flux that this function returns is SI: J s^{-1} m^{-2} m^{-1}
 *          . lambda should be between 199.5 nm and 10075 nm.
 *          . This function is overloaded.
 * 
 * @param lambda wavelength (in SI units: m)
 * @return solar radiant flux at air mass zero (in SI units: J s^{-1} m^{-2} m^{-1})
 */
double Sky::SolarRadiantFlux ( const double lambda )
{

   int i;         // Index of the first point, defining the linear relation
//   int m;         // Index of middle point
   double lam;    // wavelength [nm]
   double result; // solar radiant flux in SI units.

   // The tabulated wavelengths are in nm, so we have to convert
   // from [m] to [nm]

   lam = lambda * 1.0e+9;

   // Find the location of the wavelength point lam in the wavelength
   // array skydata::sunwavel[]

   locate ( lam, skydata::sunwavel, 920, i );

   // Check if the location was found

   if ( i == -1 )
   {
      std::cerr << "\nError (Sky::SolarRadiantFlux()): no data in this part of "
      << "the spectrum" << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }

   // Do the linear interpolation

   result = skydata::sunflux[i]
            + ( skydata::sunflux[i+1] - skydata::sunflux[i] )
            / ( skydata::sunwavel[i+1] - skydata::sunwavel[i] )
            * ( lam - skydata::sunwavel[i] );

   // The tabulated flux data are per nm (unit of wavelength), so we should
   // convert to flux per m.

   return ( result * 1.0e9 );             // conversion: (nm)^{-1} -> m^{-1}


} // end SolarRadiantFlux()
//==============================================================================




//==============================================================================
/**
 * Computes the solar radiant flux between the wavelengths lambda1 and lambda2.
 * REMARKS: . This function uses the Wehrli (1985) solar irradiance table
 *            plus the SolarRadiantFlux(lambda) function.
 *          . lambda1 and lambda2 should be between 199.5 nm and 10075 nm.
 *          . This function is overloaded.
 * 
 * @param lambda1 lower wavelength [m] of the interval
 * @param lambda2 upper wavelength [m] of the interval
 * @return solar radiant flux in [lambda1, lambda2] in J s^{-1} m^{-2}

 */
double Sky::SolarRadiantFlux ( const double lambda1,
                               const double lambda2 )
{
   double lam1;
   double lam2;

   // Check if the wavelengths are within the table boundaries.

   if (   ( lambda1 >= 199.5e-9 ) && ( lambda1 <= 10075.0e-9 )
          && ( lambda2 >= 199.5e-9 ) && ( lambda2 <= 10075.0e-9 )
      )
   {
      // If the integral boundaries are equal, the integral is zero

      if ( lambda1 == lambda2 )
      {
         return ( 0.0 );
      }

      // Check if the first wavelength is indeed greater than the second

      if ( lambda1 < lambda2 )
      {
         lam1 = lambda1;
         lam2 = lambda2;
      }
      else
      {
         lam1 = lambda2;
         lam2 = lambda1;
      }

      // Integrate with the trapezium method (Numerical Recipes, pg. 137)

      const int JMAX = 30;
      const double EPS = 1.0e-5;
      double x, tnm, sum, del, olds;
      double s = 0.0;
      int it, j, k;

      olds = -1.0e30;

      for ( j = 1; j <= JMAX; j++ )
      {
         if ( j == 1 )
         {
            s = 0.5 * ( lam2 - lam1 )
                * ( SolarRadiantFlux ( lam1 ) + SolarRadiantFlux ( lam2 ) );
         }
         else
         {
            for ( it = 1, k = 1; k < j - 1; k++ ) it <<= 1;
            tnm = it;
            del = ( lam2 - lam1 ) / tnm;
            x = lam1 + 0.5 * del;
            for ( sum = 0.0, k = 1; k <= it; k++, x += del )
            {
               sum += SolarRadiantFlux ( x );
            }
            s = 0.5 * ( s + ( lam2 - lam1 ) * sum / tnm );
         }

         if ( j > 5 )
         {
            if ( fabs ( s - olds ) < EPS * fabs ( olds ) ||
                  ( s == 0.0 && olds == 0.0 ) ) 
                return s;
         }

         olds = s;
      }

      std::cerr << "\nError (Sky::SolarRadiantFlux()): Integration not converged"<<std::endl;
      exit (1);
   }
   else
   {
      // This is the case that the given wavelengths were not between
      // the table boundaries.

      std::cerr << "\nError (Sky::SolarRadiantFlux()): wavelength must be in "
           << "[199.5e-9, 10075.0e-9]"<<std::endl;
      exit (1);
   }

   exit(0);

} // end SolarRadiantFlux()
//==============================================================================




//==============================================================================
/**
 * Compute the Solar radiant flux in the given passband
 * REMARK: . This function uses the Wehrli (1985) solar irradiance table
 *            plus the SolarRadiantFlux(lambda) function.
 *         . the passband wavelengths should be between 199.5 nm and 10075 nm.
 *         . This function is overloaded. 
 * 
 * @param lambda wavelengths of the passband [m]
 * 
 * @param throughput throughput of the passband
 * @return solar radiant flux   [J s^{-1} m^{-2}]
 */
double Sky::SolarRadiantFlux ( vector<double> &lambda,
                               vector<double> &throughput )
{
   double lambda1 = lambda[0];
   double lambda2 = lambda[lambda.size()-1];
   double integral;
//        . alpha  : right ascension coordinate  [radians]
//        . delta  : declination coordinate      [radians]
   // Check if the passband wavelengths are within [199.5, 10075] nm.

   if (   ( lambda1 < 199.5e-9 ) || ( lambda1 > 10075.0e-9 )
          || ( lambda2 < 199.5e-9 ) || ( lambda2 > 10075.0e-9 ) )
   {

      std::cerr << "\nError (Sky::SolarRadiantFlux()): Passband wavelengths not in [199.5, 10075] nm." << std::endl;
      exit (1);
   }


   // Build up the function you want to integrate

   integrand.clear();
   integrand.resize ( lambda.size() );

   for ( unsigned int i = 0; i < lambda.size(); i++ )
   {
      integrand[i] = SolarRadiantFlux ( lambda[i] ) * throughput[i];
   }

   tabfunction.Init ( lambda, integrand, lambda.size() );
   tabfunction.Set_Interpolation_Method ( Linear_Interpolation );

   integral = tabfunction.Integrate ( lambda[0], lambda[lambda.size()-1] );

   return ( integral );


} // SolarRadiantFlux()
//==============================================================================




//==============================================================================
/**
 * Compute the zodiacal background flux in the wavelength interval 
 * [lambda1, lambda2], for a given position in the sky.
 * REMARK: Some parts of the sky cannot be sampled!
 * 
 * @param lambda1 lower wavelength of the interval [m]
 * @param lambda2 upper wavelength of the interval [m]
 * @param alpha right ascension coordinate  [radians]
 * @param delta declination coordinate      [radians]
 * @return zodiacal flux [J s^{-1} m^{-2} sr^{-1}]
 */
double Sky::ZodiacalFlux ( const double lambda1, const double lambda2,
                           const double alpha, const double delta )
{
   double lam, beta;
   double flux500;
   int lam_index, beta_index;


   // Convert the equatorial coordinates to geocentric ecliptic coordinates.
   // All coordinates are in radians.

   equa2ecl ( alpha, delta, lam, beta );

   // The zodiacal light is approximately symmetric with respect to the
   // ecliptic, and with respect to the helioecliptic meridian
   // (= sun-ecliptic poles-antisolar point).

   beta = fabs ( beta );
   if ( lam > Constants::Pi ) lam = 2.0 * Constants::Pi - lam;

   // Convert from radians to degrees

   beta *= Constants::RAD2DEG;
   lam  *= Constants::RAD2DEG;

   // Locate the coordinates lam & beta in the coordinate arrays
   // Check if the coordinates are out of boundary.

   locate ( lam,  skydata::zodlong, 19, lam_index );
   locate ( beta, skydata::zodlat,  10, beta_index );

   if ( ( lam_index == -1 ) || ( beta_index == -1 ) )
   {
      std::cerr << "\nError (Sky::ZodiacalFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }

   // Check if we don't happen to be in a "hole" in the table

   if (   ( skydata::zod[lam_index][beta_index] == -1 )
          || ( skydata::zod[lam_index][beta_index+1] == -1 )
          || ( skydata::zod[lam_index+1][beta_index] == -1 )
          || ( skydata::zod[lam_index+1][beta_index+1] == -1 )
      )
   {
      std::cerr << "\nError (Sky::ZodiacalFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }

   // Do a bilinear interpolation

   double dx = ( skydata::zodlat[beta_index+1] - beta )
               / ( skydata::zodlat[beta_index+1] - skydata::zodlat[beta_index] );
   double dy = ( skydata::zodlong[lam_index+1] - lam )
               / ( skydata::zodlong[lam_index+1] - skydata::zodlong[lam_index] );
   double P1 = skydata::zod[lam_index][beta_index];
   double P2 = skydata::zod[lam_index][beta_index+1];
   double P3 = skydata::zod[lam_index+1][beta_index];
   double P4 = skydata::zod[lam_index+1][beta_index+1];

   flux500 =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
               + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;

   // The tabulated fluxes are given in:
   //    10^{-8} J s^{-1} m^{-2} sr^{-1} (micrometer)^{-1}
   // so we need to convert to:
   //    J s^{-1} m^{-2} sr^{-1} m^{-1}

   flux500 *= 0.01;

   // Now we have the zodiacal light at 500 nm, from which we need to
   // derive the zodiacal light flux in the interval [lambda1, lambda2].
   // For this we use the fact that the zodiacal flux has a solar
   // wavelength dependence.

   return ( flux500 * SolarRadiantFlux ( lambda1, lambda2 )
            / SolarRadiantFlux ( 500e-9 ) );


} // end ZodiacalFlux()
//==============================================================================




//==============================================================================
/**
 * Compute the zodiacal background flux in the given passband, 
 * for a given position in the sky.
 * REMARK: Some parts of the sky cannot be sampled!
 * 
 * @param lambda wavelengths of the passband [m]
 * @param throughput throughput of the passband
 * @param alpha right ascension coordinate  [radians]
 * @param delta declination coordinate      [radians]
 * @return  zodiacal flux [J s^{-1} m^{-2} sr^{-1}]
 */
double Sky::ZodiacalFlux ( vector<double> &lambda,
                           vector<double> &throughput,
                           const double alpha, const double delta )
{
   double lam, beta;
   double flux500;
   int lam_index, beta_index;


   // Convert the equatorial coordinates to geocentric ecliptic coordinates.
   // All coordinates are in radians.

   equa2ecl ( alpha, delta, lam, beta );

   // The zodiacal light is approximately symmetric with respect to the
   // ecliptic, and with respect to the helioecliptic meridian
   // (= sun-ecliptic poles-antisolar point).

   beta = fabs ( beta );
   if ( lam > Constants::Pi ) lam = 2.0 * Constants::Pi - lam;

   // Convert from radians to degrees

   beta *= Constants::RAD2DEG;
   lam  *= Constants::RAD2DEG;

   // Locate the coordinates lam & beta in the coordinate arrays
   // Check if the coordinates are out of boundary.

   locate ( lam,  skydata::zodlong, 19, lam_index );
   locate ( beta, skydata::zodlat,  10, beta_index );

   if ( ( lam_index == -1 ) || ( beta_index == -1 ) )
   {
      std::cerr << "\nError (Sky::ZodiacalFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }

   // Check if we don't happen to be in a "hole" in the table

   if (   ( skydata::zod[lam_index][beta_index] == -1 )
          || ( skydata::zod[lam_index][beta_index+1] == -1 )
          || ( skydata::zod[lam_index+1][beta_index] == -1 )
          || ( skydata::zod[lam_index+1][beta_index+1] == -1 )
      )
   {
      std::cerr << "\nError (Sky::ZodiacalFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }

   // Do a bilinear interpolation

   double dx = ( skydata::zodlat[beta_index+1] - beta )
               / ( skydata::zodlat[beta_index+1] - skydata::zodlat[beta_index] );
   double dy = ( skydata::zodlong[lam_index+1] - lam )
               / ( skydata::zodlong[lam_index+1] - skydata::zodlong[lam_index] );
   double P1 = skydata::zod[lam_index][beta_index];
   double P2 = skydata::zod[lam_index][beta_index+1];
   double P3 = skydata::zod[lam_index+1][beta_index];
   double P4 = skydata::zod[lam_index+1][beta_index+1];

   flux500 =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
               + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;

   // The tabulated fluxes are given in:
   //    10^{-8} J s^{-1} m^{-2} sr^{-1} (micrometer)^{-1}
   // so we need to convert to:
   //    J s^{-1} m^{-2} sr^{-1} m^{-1}

   flux500 *= 0.01;

   // Now we have the zodiacal light at 500 nm, from which we need to
   // derive the zodiacal light flux in the passband 'throughput'.
   // For this we use the fact that the zodiacal flux has a solar
   // wavelength dependence.

   return ( flux500 * SolarRadiantFlux ( lambda, throughput )
            / SolarRadiantFlux ( 500e-9 ) );



} // ZodiacalFlux()
//==============================================================================




//==============================================================================
/**
 * Return the Stellar background (unresolved stars + diffuse 
 * galactic background + extragalactic background) for the 
 * given equatorial coordinates, in the wavelength interval [lambda1, lambda2].
 * REMARK: . Some parts of the sky cannot be sampled!
 *         . We have in fact only tabulated values for the Pioneer 10
 *           blue passband (spectral range: [395, 495] nm), and for the
 *           Pioneer 10 red passband (spectral range: [590, 690] nm).
 *           The value in the given interval will be computer by linear
 *           inter- and extrapolation which, I hope, should not be totally
 *           wrong in the visual spectrum. Care is taken that in the Pioneer
 *           10 passbands, the interpolation returns exactly the tabulated
 *           values. 
 * 
 * @param lambda1 begin wavelength of the interval [m]
 * @param lambda2 end   wavelength of the interval [m]
 * @param RA equatorial coordinate: right ascension [radians]
 * @param decl equatorial coordinate: declination [radians]
 * @return Stellar background flux in the Pioneer 10 blue/red passband [J s^{-1} m^{-2} sr^{-1}]
 */
double Sky::StellarBgFlux ( const double lambda1, const double lambda2,
                            const double RA, const double decl )
{
   double alpha, delta;
   int alpha_index, delta_index;
   double blueflux, redflux;
   double a, b;

   // Convert from radians to degrees

   alpha = RA * Constants::RAD2DEG;
   delta = decl * Constants::RAD2DEG;

   // Locate the coordinates alpha & delta in the coordinate arrays
   // Check if the coordinates are out of boundary.

   locate ( alpha, skydata::skyRA, 37, alpha_index );
   locate ( delta, skydata::skydec, 25, delta_index );

   if ( ( alpha_index == -1 ) || ( delta_index == -1 ) )
   {
      std::cerr << "\nError (Sky::StellarBgFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }


   // Check if we don't happen to be in a "hole" in the table
   // These "holes" are the same for the blue and red passband.

   if (   ( skydata::skyblue[alpha_index][delta_index] == -1 )
          || ( skydata::skyblue[alpha_index][delta_index+1] == -1 )
          || ( skydata::skyblue[alpha_index+1][delta_index] == -1 )
          || ( skydata::skyblue[alpha_index+1][delta_index] == -1 )
      )
   {
      std::cerr << "\nError (Sky::StellarBgFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }


   // Do a bilinear interpolation for both the blue and the red passband

   double dx = ( skydata::skydec[delta_index+1] - delta )
               / ( skydata::skydec[delta_index+1] - skydata::skydec[delta_index] );
   double dy = ( skydata::skyRA[alpha_index+1] - alpha )
               / ( skydata::skyRA[alpha_index+1] - skydata::skyRA[alpha_index] );

   double P1 = skydata::skyblue[alpha_index][delta_index];
   double P2 = skydata::skyblue[alpha_index][delta_index+1];
   double P3 = skydata::skyblue[alpha_index+1][delta_index];
   double P4 = skydata::skyblue[alpha_index+1][delta_index+1];

   blueflux =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
                + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;

   P1 = skydata::skyred[alpha_index][delta_index];
   P2 = skydata::skyred[alpha_index][delta_index+1];
   P3 = skydata::skyred[alpha_index+1][delta_index];
   P4 = skydata::skyred[alpha_index+1][delta_index+1];

   redflux =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
               + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;


   // Convert the result from S10sun units (brightness of 10th magnitude
   // solar type stars per degree square) to SI units:
   //   J s^{-1} m^{-2} sr^{-1}

   blueflux *= 1.2084e-9;
   redflux  *= 1.0757e-9;


   // Do an extra/interpolation to the user given wavelength band
   // Note: Our 'a' is in the formulae 'a/2'.


   const double R1 = 590.0e-9;   // Lower wavelength of Pioneer 10 Red Passb.
   const double R2 = 690.0e-9;   // Upper wavelength of Pioneer 10 Red Passb.
   const double B1 = 395.0e-9;   // Lower wavelength of Pioneer 10 Blue Passb.
   const double B2 = 495.0e-9;   // Upper wavelength of Pioneer 10 Blue Passb.

   a =  ( redflux * ( B2 - B1 ) - blueflux * ( R2 - R1 ) )
        / ( ( R2 * R2 - R1 * R1 ) * ( B2 - B1 ) - ( B2 * B2 - B1 * B1 ) * ( R2 - R1 ) );

   b =  ( blueflux * ( R2 * R2 - R1 * R1 ) - redflux * ( B2 * B2 - B1 * B1 ) )
        / ( ( R2 * R2 - R1 * R1 ) * ( B2 - B1 ) - ( B2 * B2 - B1 * B1 ) * ( R2 - R1 ) );


   return (  a * ( lambda2*lambda2 - lambda1*lambda1 )
             + b * ( lambda2 - lambda1 ) );


} // end StellarBgFlux()
//==============================================================================




//==============================================================================
/**
 * Return a rough approximation of the stellar background 
 * (unresolved stars + diffuse galactic background + extragalactic background) 
 * for the given equatorial coordinates, in the given passband.
 * REMARK: . Some parts of the sky cannot be sampled!
 *         . We have in fact only tabulated values for the Pioneer 10
 *           blue passband (spectral range: [395, 495] nm), and for the
 *           Pioneer 10 red passband (spectral range: [590, 690] nm).
 *         . For lambda <= 690 nm, first the monochromatic flux function
 *           is approximated as a linear function (increasing with wavelength
 *           because there is more red than blue background light). Care is
 *           taken that the approximated flux is never negative. For
 *           lambda >= 690 nm, the monochromatic flux is approximated as
 *           a constant function with the value of the flux at lambda = 690 nm.
 *         . Care is taken that in the Pioneer 10 passbands, the interpolation
 *           returns exactly the tabulated values.
 *   
 * @param lambda wavelength values of the passband  [m]
 * @param throughput throughput of the passband
 * @param RA equatorial coordinate: right ascension [radians]
 * @param decl equatorial coordinate: declination [radians]
 * @return  A rough estimate of the stellar background flux in the given 
 * passband, given the flux in the Pioneer 10 blue/red passbands. [J s^{-1} m^{-2} sr^{-1}]
 */
double Sky::StellarBgFlux ( vector<double> &lambda,
                            vector<double> &throughput,
                            const double RA, const double decl )
{
   double alpha, delta;
   int alpha_index, delta_index;
   double blueflux, redflux;
   double a, b;

   // Convert from radians to degrees

   alpha = RA * Constants::RAD2DEG;
   delta = decl * Constants::RAD2DEG;

   // Locate the coordinates alpha & delta in the coordinate arrays
   // Check if the coordinates are out of boundary.

   locate ( alpha, skydata::skyRA, 37, alpha_index );
   locate ( delta, skydata::skydec, 25, delta_index );

   if ( ( alpha_index == -1 ) || ( delta_index == -1 ) )
   {
      std::cerr << "\nError (Sky::StellarBgFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }


   // Check if we don't happen to be in a "hole" in the table
   // These "holes" are the same for the blue and red passband.

   if (   ( skydata::skyblue[alpha_index][delta_index] == -1 )
          || ( skydata::skyblue[alpha_index][delta_index+1] == -1 )
          || ( skydata::skyblue[alpha_index+1][delta_index] == -1 )
          || ( skydata::skyblue[alpha_index+1][delta_index] == -1 )
      )
   {
      std::cerr << "\nError (Sky::StellarBgFlux()): No data for this part of the sky." << std::endl;
      std::cerr << "Please enter a value for the background flux manually in the parameter file." << std::endl;
      exit ( 1 );
   }


   // Do a bilinear interpolation for both the blue and the red passband

   double dx = ( skydata::skydec[delta_index+1] - delta )
               / ( skydata::skydec[delta_index+1] - skydata::skydec[delta_index] );
   double dy = ( skydata::skyRA[alpha_index+1] - alpha )
               / ( skydata::skyRA[alpha_index+1] - skydata::skyRA[alpha_index] );

   double P1 = skydata::skyblue[alpha_index][delta_index];
   double P2 = skydata::skyblue[alpha_index][delta_index+1];
   double P3 = skydata::skyblue[alpha_index+1][delta_index];
   double P4 = skydata::skyblue[alpha_index+1][delta_index+1];

   blueflux =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
                + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;

   P1 = skydata::skyred[alpha_index][delta_index];
   P2 = skydata::skyred[alpha_index][delta_index+1];
   P3 = skydata::skyred[alpha_index+1][delta_index];
   P4 = skydata::skyred[alpha_index+1][delta_index+1];

   redflux =   dx * dy * P1 + ( 1.0 - dx ) * dy * P2
               + dx * ( 1.0 - dy ) * P3  + ( 1.0 - dx ) * ( 1.0 - dy ) * P4;


   // Convert the result from S10sun units (brightness of 10th magnitude
   // solar type stars per degree square) to SI units:
   //   J s^{-1} m^{-2} sr^{-1}

   blueflux *= 1.2084e-9;
   redflux  *= 1.0757e-9;

   // Compute the linear relation f(lambda) = a * lambda + b
   // to be used for lambda <= 690 nm.

   const double R1 = 590.0e-9;   // Lower wavelength of Pioneer 10 Red Passb.
   const double R2 = 690.0e-9;   // Upper wavelength of Pioneer 10 Red Passb.
   const double B1 = 395.0e-9;   // Lower wavelength of Pioneer 10 Blue Passb.
   const double B2 = 495.0e-9;   // Upper wavelength of Pioneer 10 Blue Passb.

   a =  2.0 * ( redflux * ( B2 - B1 ) - blueflux * ( R2 - R1 ) )
        / ( ( R2 * R2 - R1 * R1 ) * ( B2 - B1 ) - ( B2 * B2 - B1 * B1 ) * ( R2 - R1 ) );

   b =  ( blueflux * ( R2 * R2 - R1 * R1 ) - redflux * ( B2 * B2 - B1 * B1 ) )
        / ( ( R2 * R2 - R1 * R1 ) * ( B2 - B1 ) - ( B2 * B2 - B1 * B1 ) * ( R2 - R1 ) );


   // Check if the monochromatic flux is indeed increasing for
   // lambda <= 690 nm.

   if ( a <= 0.0 )
   {
       LogManager::log << "WARNING: Sky::StellarBgFlux(): Unexpected behaviour of "
                       << "monochromatic background flux for lambda <= 690 nm." << std::endl;
       GlobalVariables::logManager.LogManagerShowLog(); 

   }

   // Compute where this linear function is zero. Useful to know, because
   // we don't want negative sky background fluxes.

   double root = - b / a;

   // Set up the integrand

   integrand.clear();
   integrand.resize ( lambda.size() );

   for ( unsigned int i = 0; i < lambda.size(); i++ )
   {
      if ( lambda[i] <= root )
      {
         integrand[i] = 0.0;
      }

      if ( ( lambda[i] > root ) && ( lambda[i] <= R2 ) )
      {
         integrand[i] = ( a * lambda[i] + b ) * throughput[i];
      }

      if ( lambda[i] > R2 )
      {
         integrand[i] = ( a * R2 + b ) * throughput[i];
      }
   }

   // Return the integral

   tabfunction.Init ( lambda, integrand, lambda.size() );
   tabfunction.Set_Interpolation_Method ( Linear_Interpolation );
   return ( tabfunction.Integrate ( lambda[0], lambda[lambda.size()-1] ) );


} // end StellarBgFlux()
//==============================================================================




//==============================================================================
/**
 * Given an array[0..N-1] of ascending values, and a value x, this 
 * function returns an index so that 
 * array[index] <= x <= array[index+1]
 * If the value x is out of the array boundaries, index will be set to -1.
 * REMARK: . NO error trapping!
 * 
 * @param x
 * @param array
 * @param N
 * @param index
 */
void Sky::locate ( double x, const double *array, int N, int &index )
{
   int index1, index2;

   // Check if x is out of the tabulated range

   if ( ( x < array[0] ) || ( x > array[N-1] ) )
   {
      index = -1;
      return;
   }

   // Check if x happens to be the first element of the array

   if ( x == array[0] )
   {
      index = 0;
   }

   // Check if x happens to be the last element of the array

   if ( x == array[N-1] )
   {
      index = N - 2;
   }

   // Find the location with bisection

   index1 = 0;                       // We already checked the lower and upper
   index2 = N - 1;                   // borders.

   unsigned int middle;               // Middle point

   while ( index2 - index1 > 1 )
   {
      middle = ( index2 + index1 ) >> 1;

      if ( x >= array[middle] )
      {
         index1 = middle;
      }
      else
      {
         index2 = middle;
      }
   }

   index = index1;

} // end locate()
//==============================================================================




//==============================================================================
/**
 * convert from equatorial coordinates (alpha = right ascension,
 *  delta = declination) to geocentric ecliptic coordinates (lambda = longitude, beta = latitude)
 * REMARK: for formulas see: Leinert et al., 1998, A&ASS 127, pg. 5
 *         for figure   see: Bernstein et al., 2002, ApJ 571, pg. 86
 * @param alpha right ascension [radians]
 * @param delta declination     [radians]
 * @param lambda ecliptic longitude [radians]
 * @param beta ecliptic latitude  [radians]
 */
void Sky::equa2ecl ( const double alpha, const double delta,
                     double &lambda, double &beta )
{
   double sinbeta, cosbeta;
   double coslambda, sinlambda;


   sinbeta = sin ( delta ) * cos ( obliquity )
             - cos ( delta ) * sin ( obliquity ) * sin ( alpha );

   beta = asin ( sinbeta );
   cosbeta = cos ( beta );

   if ( cosbeta == 0.0 )
   {
      std::cerr << "\nError (Sky::equa2ecl()): pointing ecliptic pole." << std::endl;
      exit (1);
   }

   sinlambda = ( sin ( delta ) * sin ( obliquity )
                 + cos ( delta ) * cos ( obliquity ) * sin ( alpha ) ) / cosbeta;

   coslambda = cos ( alpha ) * cos ( delta ) / cosbeta;

   lambda = atan2 ( sinlambda, coslambda );

   if ( lambda < 0.0 ) lambda += 2.0 * Constants::Pi;

} // end equa2ecl()
//==============================================================================


