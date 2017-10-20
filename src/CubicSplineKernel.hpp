/*******************************************************************************
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
 * @file CubicSplineKernel.hpp
 *
 * @brief Cubic spline kernel function.
 *
 * Based on the expression implemented in Gadget2 (Springel, 2005).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CUBICSPLINEKERNEL_HPP
#define CUBICSPLINEKERNEL_HPP

/**
 * @brief Cubic spline kernel function.
 *
 * Based on the expression implemented in Gadget2 (Springel, 2005).
 */
class CubicSplineKernel {
public:
  /**
   * @brief Evaluate the kernel at the given normalized distance.
   *
   * @param u Normalized distance (in units of h).
   * @param h Smoothing length (in m).
   * @return Value of the cubic spline kernel.
   */
  static inline double kernel_evaluate(double u, double h) {

    const double KC1 = 2.546479089470;
    const double KC2 = 15.278874536822;
    const double KC5 = 5.092958178941;
    if (u < 1.) {
      if (u < 0.5) {
        return (KC1 + KC2 * (u - 1.) * u * u) / (h * h * h);
      } else {
        return KC5 * (1. - u) * (1. - u) * (1. - u) / (h * h * h);
      }
    } else {
      // the cubic spline kernel has compact support
      return 0.;
    }
  }
};

#endif // CUBICSPLINEKERNEL_HPP
