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
 * @file IsotropicContinuousPhotonSource.hpp
 *
 * @brief Class used to generate an isotropic external radiation field.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ISOTROPICCONTINUOUSPHOTONSOURCE_HPP
#define ISOTROPICCONTINUOUSPHOTONSOURCE_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"

#include <cfloat>
#include <limits>

/**
 * @brief Class used to generate an isotropic external radiation field.
 */
class IsotropicContinuousPhotonSource {
private:
  /*! @brief Box in which the radiation enters. */
  Box _box;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box in which the radiation enters (in m).
   * @param log Log to write logging info to.
   */
  IsotropicContinuousPhotonSource(Box box, Log *log = nullptr) : _box(box) {

    if (log) {
      log->write_status(
          "Constructed IsotropicContinuousPhotonSource in box with anchor [",
          _box.get_anchor().x(), " m, ", _box.get_anchor().y(), " m, ",
          _box.get_anchor().z(), " m], and sides [", _box.get_sides().x(),
          " m, ", _box.get_sides().y(), " m, ", _box.get_sides().z(), " m].");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParamterFile to read from.
   * @param log Log to write logging info to.
   */
  IsotropicContinuousPhotonSource(ParameterFile &params, Log *log = nullptr)
      : IsotropicContinuousPhotonSource(
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides")),
            log) {}

  /**
   * @brief Get the entrance position and direction of a random external photon.
   *
   * @param random_generator RandomGenerator to use.
   * @return std::pair of CoordinateVector instances, specifying a starting
   * position and direction for an incoming photon.
   */
  std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const {
    // we randomly sample a focus point in the box
    CoordinateVector<> focus;
    focus[0] =
        _box.get_anchor().x() +
        _box.get_sides().x() * random_generator.get_uniform_random_double();
    focus[1] =
        _box.get_anchor().y() +
        _box.get_sides().y() * random_generator.get_uniform_random_double();
    focus[2] =
        _box.get_anchor().z() +
        _box.get_sides().z() * random_generator.get_uniform_random_double();

    // random incoming direction for the focus point
    double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    double sint = 1. - cost * cost;
    sint = std::sqrt(std::max(sint, 0.));
    double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    double cosp = std::cos(phi);
    double sinp = std::sin(phi);
    CoordinateVector<> direction(sint * cosp, sint * sinp, cost);

    // the direction and the focus point define a line:
    //   line = focus + t * direction.
    // with t a parameter.
    // find the intersection points of this line with the walls of the box
    // the origin of the random photon is the intersection point with negative
    // t value
    CoordinateVector<> anchor_bottom = _box.get_anchor();
    CoordinateVector<> anchor_top = _box.get_top_anchor();

    // the box has 6 faces, each of which has an intersection point with the
    // line. If the focus point is inside the box (which it should be), then
    // 3 intersection points will have a negative parameter value, and 3 will
    // have a positive value
    // we record all 3 negative values, and then pick the one with the largest
    // parameter value (smallest absolute value), since this intersection point
    // is the closest
    double lx, ly, lz;
    if (direction.x() < 0.) {
      lx = (anchor_top.x() - focus.x()) / direction.x();
    } else if (direction.x() > 0.) {
      lx = (anchor_bottom.x() - focus.x()) / direction.x();
    } else {
      // the line is parallel to the x faces
      // make sure this intersection point has the largest absolute parameter
      // value
      lx = -DBL_MAX;
    }

    if (direction.y() < 0.) {
      ly = (anchor_top.y() - focus.y()) / direction.y();
    } else if (direction.y() > 0.) {
      ly = (anchor_bottom.y() - focus.y()) / direction.y();
    } else {
      ly = -DBL_MAX;
    }

    if (direction.z() < 0.) {
      lz = (anchor_top.z() - focus.z()) / direction.z();
    } else if (direction.z() > 0.) {
      lz = (anchor_bottom.z() - focus.z()) / direction.z();
    } else {
      lz = -DBL_MAX;
    }

    double maxl = std::max(lx, ly);
    maxl = std::max(maxl, lz);

    CoordinateVector<> position = focus + maxl * direction;

    // make sure the photon is inside the box
    for (unsigned int i = 0; i < 3; ++i) {
      // we cannot simply take the top anchor of the box as upper limit, since
      // the top anchor itself strictly speaking lies outside the box (lower
      // limits are inclusive, upper limits exclusive due to the way we
      // calculate grid indices)
      // we therefore take the closest value that is still in the box
      // epsilon is the difference between 1.0 and the next floating point value
      // larger than 1.0 that can be represented as a 64-bit floating point.
      position[i] = std::min(position[i],
                             anchor_top[i] -
                                 std::numeric_limits< double >::epsilon() *
                                     _box.get_sides()[i]);
      position[i] = std::max(position[i], anchor_bottom[i]);
    }

    return std::make_pair(position, direction);
  }

  /**
   * @brief Get the total surface area through which the radiation enters the
   * simulation box.
   *
   * @return Total surface area (in m^2).
   */
  inline double get_total_surface_area() const {
    return 2. * _box.get_sides().x() * _box.get_sides().y() +
           2. * _box.get_sides().x() * _box.get_sides().z() +
           2. * _box.get_sides().y() * _box.get_sides().z();
  }
};

#endif // ISOTROPICCONTINUOUSPHOTONSOURCE_HPP
