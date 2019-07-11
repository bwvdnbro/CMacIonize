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
 * @file RateBasedUVLuminosityFunction.hpp
 *
 * @brief UVLuminosityFunction implementation that uses a fixed UV rate per mass
 * unit.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RATEBASEDUVLUMINOSITYFUNCTION_HPP
#define RATEBASEDUVLUMINOSITYFUNCTION_HPP

#include "ParameterFile.hpp"
#include "UVLuminosityFunction.hpp"

/**
 * @brief UVLuminosityFunction implementation that uses a fixed UV rate per mass
 * unit.
 */
class RateBasedUVLuminosityFunction : public UVLuminosityFunction {
private:
  /*! @brief UV luminosity per unit mass (in s^-1 kg^-1). */
  const double _UV_rate_per_mass_unit;

  /*! @brief Upper limit on the age (in s). */
  const double _cutoff_age;

public:
  /**
   * @brief Constructor.
   *
   * @param UV_rate_per_mass_unit UV luminosity per mass unit (in s^-1 kg^-1).
   * @param cutoff_age Upper limit on the age (in s).
   */
  inline RateBasedUVLuminosityFunction(const double UV_rate_per_mass_unit,
                                       const double cutoff_age)
      : _UV_rate_per_mass_unit(UV_rate_per_mass_unit), _cutoff_age(cutoff_age) {
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - UV rate per mass unit: UV luminosity per unit mass (default:
   *    2.49428e16 s^-1 kg^-1)
   *  - cutoff age: Upper limit on the age (default: 5. Myr)
   *
   * @param params ParameterFile to read from.
   */
  inline RateBasedUVLuminosityFunction(ParameterFile &params)
      : RateBasedUVLuminosityFunction(
            params.get_physical_value< QUANTITY_FREQUENCY_PER_MASS >(
                "UVLuminosityFunction:UV rate per mass unit",
                "2.49428e16 s^-1 kg^-1"),
            params.get_physical_value< QUANTITY_TIME >(
                "UVLuminosityFunction:cutoff age", "5. Myr")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~RateBasedUVLuminosityFunction() {}

  /**
   * @brief Function that gives the UV luminosity for a star particle with the
   * given age and mass.
   *
   * @param age Age of the star particle (in s).
   * @param mass Total mass of the star particle (in kg).
   * @return UV luminosity of the star particle (in s^-1).
   */
  virtual double operator()(const double age, const double mass) const {
    if (age <= _cutoff_age) {
      return mass * _UV_rate_per_mass_unit;
    } else {
      return 0.;
    }
  }
};

#endif // RATEBASEDUVLUMINOSITYFUNCTION_HPP
