/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DiscICDensityFunction.hpp
 *
 * @brief Disc initial condition DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCICDENSITYFUNCTION_HPP
#define DISCICDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction that returns a constant density value for all
 * coordinates, corresponding to a homogeneous density field, plus a cylindrical
 * tangential \f$\frac{1}{r}\f$ velocity profile.
 */
class DiscICDensityFunction : public DensityFunction {
private:
  /*! @brief Single density value for the entire box (in m^-3). */
  const double _density;

  /*! @brief Single temperature value for the entire box (in K). */
  const double _temperature;

  /*! @brief Initial hydrogen neutral fraction for the entire box. */
  const double _neutral_fraction_H;

  /*! @brief Characteristic radius for the \f$\frac{1}{r}\f$ velocity profile
   *  (in m). */
  const double _r_C;

  /*! @brief Characteristic velocity for the \f$\frac{1}{r}\f$ velocity
   *  profile (in m s^-1). */
  const double _v_C;

public:
  /**
   * @brief Constructor.
   *
   * @param density Single density value for the entire box (in m^-3).
   * @param temperature Single temperature value for the entire box (in K).
   * @param neutral_fraction_H Single hydrogen neutral fraction value for the
   * entire box.
   * @param r_C Characteristic radius for the velocity profile (in m).
   * @param v_C Characteristic velocity for the velocity profile (in m s^-1).
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(double density = 1., double temperature = 8000.,
                        double neutral_fraction_H = 1.e-6, double r_C = 1.,
                        double v_C = 1., Log *log = nullptr)
      : _density(density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H), _r_C(r_C), _v_C(v_C) {

    if (log) {
      log->write_status(
          "Created HomogeneousDensityFunction with constant density ", _density,
          " m^-3 and constant temperature ", _temperature, " K.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - density: Constant number density value (default: 100. cm^-3)
   *  - temperature: Constant initial temperature value (default: 8000. K)
   *  - neutral fraction H: Contant initial neutral fraction value
   *    (default: 1.e-6)
   *  - characteristic radius: radial parameter for the velocity profile
   *    (default: 1. m)
   *  - characteristic velocity: velocity parameter for the velocity profile
   *    (default: 1. m s^-1)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(ParameterFile &params, Log *log = nullptr)
      : DiscICDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:density", "100. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "8000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H",
                                       1.e-6),
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:characteristic radius", "1. m"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "DensityFunction:characteristic velocity", "1. m s^-1"),
            log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {

    DensityValues values;

    values.set_number_density(_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction_H);
    values.set_ionic_fraction(ION_He_n, 1.e-6);

    /// velocity field
    // get the cell position
    const CoordinateVector<> p = cell.get_cell_midpoint();
    // get the inverse cylindrical radius
    const double Rinv2 = 1. / (p.x() * p.x() + p.y() * p.y());
    // get the velocity
    const double vphi = _v_C * _r_C * Rinv2;

    const CoordinateVector<> velocity(-vphi * p.y(), vphi * p.x(), 0.);
    values.set_velocity(velocity);

    return values;
  }
};

#endif // HOMOGENEOUSDENSITYFUNCTION_HPP
