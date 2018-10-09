/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file CaproniStellarRoutines.hpp
 *
 * @brief Common routines used in CaproniPhotonSourceDistribution and
 * CaproniStellarFeedback.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CAPRONISTELLARROUTINES_HPP
#define CAPRONISTELLARROUTINES_HPP

#include "CoordinateVector.hpp"
#include "RandomGenerator.hpp"

#include <cmath>

/**
 * @brief Common routines used in CaproniPhotonSourceDistribution and
 * CaproniStellarFeedback.
 */
class CaproniStellarRoutines {
public:
  /**
   * @brief Get the expected number of stars for the given time.
   *
   * This number is based on the polynomial fit to the Caproni et al. (2017)
   * derived OB number function.
   *
   * @param t Current simulation time (in s).
   * @return Expected number of OB stars at this time.
   */
  inline static uint_fast32_t get_number_of_stars(const double t) {

    // coefficients for the polynomial
    const double a[10] = {-4.22241876763e-146, 1.66204113118e-128,
                          -2.68754635204e-111, 2.33363366028e-94,
                          -1.18287541821e-77,  3.52485869869e-61,
                          -5.7769427071e-45,   4.08847904916e-29,
                          1.01863714667e-15,   29.0068124326};

    // evaluate polynomial in O(n) time using Horner's method
    double result = a[0] * t + a[1];
    result = result * t + a[2];
    result = result * t + a[3];
    result = result * t + a[4];
    result = result * t + a[5];
    result = result * t + a[6];
    result = result * t + a[7];
    result = result * t + a[8];
    result = result * t + a[9];

    return result;
  }

  /**
   * @brief Get the UV luminosity and life time for a random OB star distributed
   * according to the IMF.
   *
   * @param luminosity Output UV luminosity of the random star (in s^-1).
   * @param lifetime Output life time of the star (in s).
   * @param random_generator RandomGenerator to use.
   * @return Mass of the star (in Msol; only for test purposes).
   */
  inline static double get_random_star(double &luminosity, double &lifetime,
                                       RandomGenerator &random_generator) {

    const double alphap1 = -1.3;
    const double alphap1inv = 1. / alphap1;
    const double mlowterm = std::pow(20., alphap1);
    const double mrangefac = std::pow(100., alphap1) - mlowterm;

    // get a random mass distributed according to the IMF
    const double u = random_generator.get_uniform_random_double();
    const double mstar = std::pow(u * mrangefac + mlowterm, alphap1inv);

    // get the corresponding luminosity
    // we use a 3th order polynomial fit (and extrapolation) for the OB
    // luminosity data from Sternberg et al. (2003) (OB_LCV.dat)
    const double a[4] = {-8.85154170718e+43, 2.21555601476e+46,
                         -4.25455875963e+47, 8.55819263554e+47};
    luminosity = a[0] * mstar + a[1];
    luminosity = luminosity * mstar + a[2];
    luminosity = luminosity * mstar + a[3];

    // get the corresponding life time
    // we use a double power law fit to the stellar life time data from
    // Tang et al. (2014), the Z0.017Y0.279 model for masses in the range
    // [20, 100] Msol. The life time was computed by taking the difference
    // between the stellar ages at the start of phase 4 (near the ZAM) and
    // phase 8 (base of the RGB).
    const double la[5] = {7.55609422e+13, 1.03371798e+16, -1.31168267e+00,
                          1.11162246e+18, -3.81030835e+00};
    lifetime =
        la[0] + la[1] * std::pow(mstar, la[2]) + la[3] * std::pow(mstar, la[4]);

    return mstar;
  }

  /**
   * @brief Get a random galactic radius corresponding to the given time.
   *
   * The radius is based on a polynomial fit to the Caproni et al. (2017) SN
   * location data.
   *
   * @param t Current simulation time (in s).
   * @param random_generator RandomGenerator to use.
   * @return Random galactic radius (in m).
   */
  inline static double get_galactic_radius(const double t,
                                           RandomGenerator &random_generator) {

    // fit coefficients;
    const double a[10] = {-2.48299482225e-128, 6.91524461551e-111,
                          -7.96019285215e-94,  4.88977912839e-77,
                          -1.72602520832e-60,  3.51127816196e-44,
                          -3.91671781513e-28,  2.10807296898e-12,
                          -2773.58198637,      4.62567627189e+18};

    double ravg = a[0] * t + a[1];
    ravg = ravg * t + a[2];
    ravg = ravg * t + a[3];
    ravg = ravg * t + a[4];
    ravg = ravg * t + a[5];
    ravg = ravg * t + a[6];
    ravg = ravg * t + a[7];
    ravg = ravg * t + a[8];
    ravg = ravg * t + a[9];

    const double gauss =
        std::sqrt(-2. *
                  std::log(random_generator.get_uniform_random_double())) *
        std::cos(2. * M_PI * random_generator.get_uniform_random_double());

    return ravg + 3.086e18 * gauss;
  }

  /**
   * @brief Generate a new source position.
   *
   * @param t Current simulation time (in s).
   * @param random_generator RandomGenerator to use.
   * @return New source position (in m).
   */
  inline static CoordinateVector<>
  generate_source_position(const double t, RandomGenerator &random_generator) {

    // get a random radius for this time
    const double r = get_galactic_radius(t, random_generator);

    // get a random direction
    const double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);

    return CoordinateVector<>(r * sint * cosp, r * sint * sinp, r * cost);
  }
};

#endif // CAPRONISTELLARROUTINES_HPP
