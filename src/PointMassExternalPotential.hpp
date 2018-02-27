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
 * @file PointMassExternalPotential.hpp
 *
 * @brief Point mass external potential.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef POINTMASSEXTERNALPOTENTIAL_HPP
#define POINTMASSEXTERNALPOTENTIAL_HPP

#include "CoordinateVector.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

/**
 * @brief Point mass external potential.
 */
class PointMassExternalPotential {
private:
  /*! @brief Position of the point mass (in m). */
  const CoordinateVector<> _position;

  /*! @brief Rescaled mass of the point mass (in units of kg * G). */
  const double _m_G;

public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the point mass (in m).
   * @param mass Mass of the point mass (in kg).
   */
  inline PointMassExternalPotential(const CoordinateVector<> position,
                                    const double mass)
      : _position(position),
        _m_G(mass * PhysicalConstants::get_physical_constant(
                        PHYSICALCONSTANT_NEWTON_CONSTANT)) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters:
   *  - position: Position of the external point mass
   *    (default: [0. m, 0. m, 0. m])
   *  - mass: Mass of the external point mass
   *
   * @param params ParameterFile to read from.
   */
  inline PointMassExternalPotential(ParameterFile &params)
      : PointMassExternalPotential(
            params.get_physical_vector< QUANTITY_LENGTH >(
                "ExternalPotential:position", "[0. m, 0. m, 0. m]"),
            params.get_physical_value< QUANTITY_MASS >("ExternalPotential:mass",
                                                       "1. kg")) {}

  /**
   * @brief Get the acceleration caused by the external point mass on a mass at
   * the given position.
   *
   * @param position Position (in m).
   * @return Acceleration (in m s^-2).
   */
  inline CoordinateVector<>
  get_acceleration(const CoordinateVector<> position) const {
    const CoordinateVector<> dx = position - _position;
    const double r2 = dx.norm2();
    const double r = std::sqrt(r2);
    return dx * (-_m_G / (r * r2));
  }
};

#endif // POINTMASSEXTERNALPOTENTIAL_HPP
