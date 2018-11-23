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
 * @file CoredDMProfileExternalPotential.hpp
 *
 * @brief Cored DM profile external potential.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COREDDMPROFILEEXTERNALPOTENTIAL_HPP
#define COREDDMPROFILEEXTERNALPOTENTIAL_HPP

#include "CoordinateVector.hpp"
#include "ExternalPotential.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Cored DM profile external potential.
 *
 * The potential is based on Caproni et al. (2015).
 */
class CoredDMProfileExternalPotential : public ExternalPotential {
private:
  /*! @brief Inverse core radius @f$\frac{1}{r_0}@f$ (in m^-1). */
  const double _r0inv;

  /*! @brief Velocity at infinity squared @f$v_\inf{}^2@f$ (in m^2 s^-2). */
  const double _vinf2;

public:
  /**
   * @brief Constructor.
   *
   * @param r0 Core radius (in m).
   * @param vinf Velocity at infinity (in m s^-1).
   */
  inline CoredDMProfileExternalPotential(const double r0, const double vinf)
      : _r0inv(1. / r0), _vinf2(vinf * vinf) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We read the following parameters:
   *  - core radius: Radius of the core of the DM profile (default: 300. pc)
   *  - maximum circular velocity: Maximum circular velocity of the DM profile
   *    (default: 21.1 km s^-1).
   *
   * @param params ParameterFile to read from.
   */
  inline CoredDMProfileExternalPotential(ParameterFile &params)
      : CoredDMProfileExternalPotential(
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:core radius", "300. pc"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "ExternalPotential:maximum circular velocity",
                "21.1 km s^-1")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~CoredDMProfileExternalPotential() {}

  /**
   * @brief Get the acceleration caused by the external disc on a mass at the
   * given position.
   *
   * @param position Position (in m).
   * @return Acceleration (in m s^-2).
   */
  virtual CoordinateVector<>
  get_acceleration(const CoordinateVector<> position) const {

    const double r = position.norm();
    const double rinv = 1. / r;
    const double ksi = r * _r0inv;
    const double ksiinv = 1. / ksi;

    const double dphidksi =
        -_vinf2 * (ksiinv - ksiinv * ksiinv * std::atan(ksi));

    return position * _r0inv * rinv * dphidksi;
  }
};

#endif // COREDDMPROFILEEXTERNALPOTENTIAL_HPP
