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
 * @file HomogeneousDensityFunction.hpp
 *
 * @brief Homogeneous DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HOMOGENEOUSDENSITYFUNCTION_HPP
#define HOMOGENEOUSDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction that returns a constant value for all coordinates,
 * corresponding to a homogeneous density field.
 */
class HomogeneousDensityFunction : public DensityFunction {
private:
  /*! @brief Single density value for the entire box. */
  double _density;

public:
  /**
   * @brief Constructor.
   *
   * @param density Single density value for the entire box.
   */
  HomogeneousDensityFunction(double density = 1.) : _density(density) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  HomogeneousDensityFunction(ParameterFile &params)
      : HomogeneousDensityFunction(
            params.get_value< double >("densityfunction.density")) {}

  /**
   * @brief Function that gives the density for a given coordinate.
   *
   * @param position CoordinateVector specifying a coordinate position.
   * @return Density at the given coordinate: the single value stored
   * internally.
   */
  virtual double operator()(CoordinateVector<> position) { return _density; }
};

#endif // HOMOGENEOUSDENSITYFUNCTION_HPP
