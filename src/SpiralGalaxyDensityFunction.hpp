/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SpiralGalaxyDensityFunction.hpp
 *
 * @brief DensityFunction implementation that returns a spiral galaxy density
 * field.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPIRALGALAXYDENSITYFUNCTION_HPP
#define SPIRALGALAXYDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction implementation that returns a spiral galaxy density
 * field.
 *
 * We assume a simple double exponential profile of the form
 * \f[
 *  n(x, y, z) = n_0 \exp\left(-\frac{\sqrt{x^2+y^2}}{r_{ISM}}\right)
 *               \exp\left(-\frac{|z|}{h_{ISM}}\right).
 * \f]
 */
class SpiralGalaxyDensityFunction : public DensityFunction {
private:
  /*! @brief Scale length @f$r_{ISM}@f$ of the gas and dust in the ISM
   *  (in m). */
  const double _r_ISM;

  /*! @brief Scale height @f$h_{ISM}@f$ of the gas and dust in the ISM
   *  (in m). */
  const double _h_ISM;

  /*! @brief Central density @f$\rho{}_0@f$ of the gas and dust in the ISM
   *  (in kg m^-3). */
  const double _n_0;

  /*! @brief Conversion factor from kpc to m (in m). */
  const double _kpc;

public:
  /**
   * @brief Constructor.
   *
   * @param rdust Scale length @f$r_{ISM}@f$ of the gas and dust in the ISM
   * (in m).
   * @param hdust Scale height @f$h_{ISM}@f$ of the gas and dust in the ISM
   * (in m).
   * @param n_0 Central density @f$n_0@f$ of the gas and dust in the ISM
   * (in m^-3).
   * @param log Log to write logging info to.
   */
  SpiralGalaxyDensityFunction(double rdust, double hdust, double n_0,
                              Log *log = nullptr)
      : _r_ISM(rdust), _h_ISM(hdust),
        /* we multiplied the value in Kenny's code with 1.e-3 to convert g to
         * kg */
        _n_0(1.674e-27 * n_0), _kpc(3.086e19) {

    if (log) {
      log->write_status(
          "Constructed SpiralGalaxyDensityFunction with scale length ", _r_ISM,
          " m, scale height ", _h_ISM, " m, and a central density of ", _n_0,
          " kg m^-3.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - scale length ISM: Scale length of the gas disc (default: 6. kpc)
   *  - scale height ISM: Scale height of the gas disc (default: 0.22 kpc)
   *  - central density: Central density of the gas disc (default: 1. cm^-3)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SpiralGalaxyDensityFunction(ParameterFile &params, Log *log = nullptr)
      : SpiralGalaxyDensityFunction(
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:scale length ISM", "6. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:scale height ISM", "0.22 kpc"),
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:central density", "1. cm^-3"),
            log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  inline DensityValues operator()(const Cell &cell) const {
    const CoordinateVector<> &position = cell.get_cell_midpoint();
    const double w2 = position.x() * position.x() + position.y() * position.y();
    const double w = std::sqrt(w2);
    double n = 0.;
    if (w < 15. * _kpc && std::abs(position.z()) < 15. * _kpc) {
      n = _n_0 * std::exp(-w / _r_ISM) *
          std::exp(-std::abs(position.z()) / _h_ISM);
    }
    DensityValues values;
    values.set_number_density(n);
    values.set_ionic_fraction(ION_H_n, 1.);
    values.set_ionic_fraction(ION_He_n, 0.);
    values.set_temperature(0.);
    return values;
  }
};

#endif // SPIRALGALAXYDENSITYFUNCTION_HPP
