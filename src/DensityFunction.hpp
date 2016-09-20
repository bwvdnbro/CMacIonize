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
 * @file DensityFunction.hpp
 *
 * @brief Interface for functors that can be used to fill a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYFUNCTION_HPP
#define DENSITYFUNCTION_HPP

class DensityFunction {
public:
  /**
   * @brief Function that gives the density for a given coordinate.
   *
   * @param x x coordinate.
   * @param y y coordinate.
   * @param z z coordinate.
   * @return Density at the given coordinate.
   */
  virtual double operator()(double x, double y, double z) = 0;
};

#endif // DENSITYFUNCTION_HPP
