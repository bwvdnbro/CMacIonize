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
 * @file UVLuminosityFunction.hpp
 *
 * @brief Interface for functors that can be used to assign UV luminosities to
 * star particles in a simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UVLUMINOSITYFUNCTION_HPP
#define UVLUMINOSITYFUNCTION_HPP

/**
 * @brief Interface for functors that can be used to assign UV luminosities to
 * star particles in a simulation.
 */
class UVLuminosityFunction {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~UVLuminosityFunction() {}

  /**
   * @brief Function that gives the UV luminosity for a star particle with the
   * given age and mass.
   *
   * @param age Age of the star particle (in s).
   * @param mass Total mass of the star particle (in kg).
   * @return UV luminosity of the star particle (in s^-1).
   */
  virtual double operator()(const double age, const double mass) const = 0;
};

#endif // UVLUMINOSITYFUNCTION_HPP
