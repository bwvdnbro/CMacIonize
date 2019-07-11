/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IMFBasedUVLuminosityFunction.hpp
 *
 * @brief UVLuminosityFunction implementation that uses an IMF to determine the
 * UV luminosity for a star particle.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IMFBASEDUVLUMINOSITYFUNCTION_HPP
#define IMFBASEDUVLUMINOSITYFUNCTION_HPP

#include "ParameterFile.hpp"
#include "UVLuminosityFunction.hpp"

#include <cmath>

/**
 * @brief UVLuminosityFunction implementation that uses an IMF to determine the
 * UV luminosity for a star particle.
 */
class IMFBasedUVLuminosityFunction : public UVLuminosityFunction {
private:
  /*! @brief High-mass slope of the IMF. */
  const double _slope;

  /*! @brief Lower mass limit for IMF integration (in Msol). */
  const double _lower_mass_limit_in_Msol;

  /*! @brief Upper mass limit for IMF integration (in Msol). */
  const double _upper_mass_limit_in_Msol;

  /*! @brief Fraction of the total IMF mass that is contained in between the
   *  lower and upper mass limit. */
  const double _high_mass_fraction;

  /*! @brief Normalisation factor for the IMF (in kg^-1). */
  const double _IMF_norm;

  /*! @brief Boost factor for UV luminosities. */
  const double _boost_factor;

  /**
   * @brief Integrate the power law tail of the IMF between the given mass
   * limits.
   *
   * @param M_L_in_Msol Lower mass limit (in Msol).
   * @param M_U_in_Msol Upper mass limit (in Msol).
   * @param slope Slope of the IMF.
   * @return Total (unnormalised) mass contained between the given mass limits.
   */
  inline static double IMF_mass_integral(const double M_L_in_Msol,
                                         const double M_U_in_Msol,
                                         const double slope) {
    return (std::pow(M_U_in_Msol, 2. - slope) -
            std::pow(M_L_in_Msol, 2. - slope)) /
           (2. - slope);
  }

  /**
   * @brief Integrate the power law tail of the IMF UV luminosity function
   * between the given mass limits.
   *
   * @param M_L_in_Msol Lower mass limit (in Msol).
   * @param M_U_in_Msol Upper mass limit (in Msol).
   * @param slope Slope of the IMF.
   * @return Total (unnormalised) UV luminosity contained between the given mass
   * limits.
   */
  inline static double IMF_UV_integral(const double M_L_in_Msol,
                                       const double M_U_in_Msol,
                                       const double slope) {

    const double A = -8.85154170718e+43;
    const double B = 2.21555601476e+46;
    const double C = -4.25455875963e+47;
    const double D = 8.55819263554e+47;

    const double p0 = (std::pow(M_U_in_Msol, 1. - slope) -
                       std::pow(M_L_in_Msol, 1. - slope)) /
                      (1. - slope);
    const double p1 = (std::pow(M_U_in_Msol, 2. - slope) -
                       std::pow(M_L_in_Msol, 2. - slope)) /
                      (2. - slope);
    const double p2 = (std::pow(M_U_in_Msol, 3. - slope) -
                       std::pow(M_L_in_Msol, 3. - slope)) /
                      (3. - slope);
    const double p3 = (std::pow(M_U_in_Msol, 4. - slope) -
                       std::pow(M_L_in_Msol, 4. - slope)) /
                      (4. - slope);

    return A * p3 + B * p2 + C * p1 + D * p0;
  }

  /**
   * @brief Get the total UV luminosity for a stellar population with the given
   * mass limit for the most massive stars.
   *
   * @param upper_mass_limit_in_Msol Mass of the most massive stars in the
   * population (in Msol).
   * @param population_mass Total mass of the entire population (in kg).
   * @return UV luminosity of the entire population (in s^-1).
   */
  inline double get_UV_luminosity(const double upper_mass_limit_in_Msol,
                                  const double population_mass) const {

    if (upper_mass_limit_in_Msol <= _lower_mass_limit_in_Msol) {
      return 0.;
    }

    const double IMF_UV_factor = IMF_UV_integral(
        _lower_mass_limit_in_Msol,
        std::min(_upper_mass_limit_in_Msol, upper_mass_limit_in_Msol), _slope);

    return IMF_UV_factor * _boost_factor * population_mass * _IMF_norm;
  }

  /**
   * @brief Get the mass of the most massive star in a stellar population with
   * the given age.
   *
   * @param age Age of the stellar population (in s).
   * @return Mass of the most massive star in the population (in Msol).
   */
  inline double get_upper_mass_limit(const double age) const {

    const double age_in_Myr = UnitConverter::convert(age, "s", "Myr");

    const double la[5] = {4.47959896e+00, 1.52686581e+02, -1.04819293e+00,
                          5.51939499e+03, -4.11097721e+00};
    return la[0] + la[1] * std::pow(age_in_Myr, la[2]) +
           la[3] * std::pow(age_in_Myr, la[4]);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param slope Slope of the IMF.
   * @param lower_mass_limit Lower mass limit for UV luminous sources (in kg).
   * @param upper_mass_limit Upper mass limit for the IMF (in kg).
   * @param high_mass_fraction Mass fraction contained in the massive end of the
   * IMF.
   * @param boost_factor Multiplicative factor to artificially enhance or reduce
   * UV luminosities.
   */
  inline IMFBasedUVLuminosityFunction(const double slope,
                                      const double lower_mass_limit,
                                      const double upper_mass_limit,
                                      const double high_mass_fraction,
                                      const double boost_factor)
      : _slope(slope),
        _lower_mass_limit_in_Msol(lower_mass_limit /
                                  PhysicalConstants::get_physical_constant(
                                      PHYSICALCONSTANT_SOLAR_MASS)),
        _upper_mass_limit_in_Msol(upper_mass_limit /
                                  PhysicalConstants::get_physical_constant(
                                      PHYSICALCONSTANT_SOLAR_MASS)),
        _high_mass_fraction(high_mass_fraction),
        _IMF_norm(_high_mass_fraction /
                  (IMF_mass_integral(_lower_mass_limit_in_Msol,
                                     _upper_mass_limit_in_Msol, _slope) *
                   PhysicalConstants::get_physical_constant(
                       PHYSICALCONSTANT_SOLAR_MASS))),
        _boost_factor(boost_factor) {

    if (_slope < 0.) {
      cmac_error("Negative slope provided for IMF (minus sign is already "
                 "assumed internally)!");
    }
    if (_slope == 1. || _slope == 2. || _slope == 3. || _slope == 4.) {
      cmac_error("Value for IMF slope leads to problems with integration: %g!",
                 _slope);
    }
    if (_boost_factor <= 0.) {
      cmac_error("Boost factor should be positive!");
    }
    if (_lower_mass_limit_in_Msol >= _upper_mass_limit_in_Msol) {
      cmac_error("Upper mass limit should be higher than lower mass limit!");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the file:
   *  - slope: Slope of the IMF (default: 2.3)
   *  - lower mass limit: Lower mass limit for UV luminous sources (default: 15.
   *    Msol)
   *  - upper mass limit: Upper mass limit for stars in the IMF (default: 100.
   *    Msol)
   *  - high mass fraction: Fraction of the IMF mass that resides in high mass
   *    stars (default: 0.25)
   *  - boost factor: Additional multiplicative factor to enhance or suppress UV
   *    luminosities (default: 1.)
   *
   * @param params ParameterFile to read from.
   */
  inline IMFBasedUVLuminosityFunction(ParameterFile &params)
      : IMFBasedUVLuminosityFunction(
            params.get_value< double >("UVLuminosityFunction:slope", 2.3),
            params.get_physical_value< QUANTITY_MASS >(
                "UVLuminosityFunction:lower mass limit", "15. Msol"),
            params.get_physical_value< QUANTITY_MASS >(
                "UVLuminosityFunction:upper mass limit", "100. Msol"),
            params.get_value< double >(
                "UVLuminosityFunction:high mass fraction", 0.25),
            params.get_value< double >("UVLuminosityFunction:boost factor",
                                       1.)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~IMFBasedUVLuminosityFunction() {}

  /**
   * @brief Function that gives the UV luminosity for a star particle with the
   * given age and mass.
   *
   * @param age Age of the star particle (in s).
   * @param mass Total mass of the star particle (in kg).
   * @return UV luminosity of the star particle (in s^-1).
   */
  virtual double operator()(const double age, const double mass) const {
    return get_UV_luminosity(get_upper_mass_limit(age), mass);
  }
};

#endif // IMFBASEDUVLUMINOSITYFUNCTION_HPP
