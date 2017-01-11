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
 * @file DistantStarContinuousPhotonSource.hpp
 *
 * @brief ContinuousPhotonSource implementation for a distant star (outside the
 * simulation box).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISTANTSTARCONTINUOUSPHOTONSOURCE_HPP
#define DISTANTSTARCONTINUOUSPHOTONSOURCE_HPP

#include "Box.hpp"
#include "ContinuousPhotonSource.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief ContinuousPhotonSource implementation for a distant star (outside the
 * simulation box).
 */
class DistantStarContinuousPhotonSource : public ContinuousPhotonSource {
private:
  /*! @brief Position of the distant star (in m). */
  CoordinateVector<> _position;

  /*! @brief Flags that keep track of which sides of the box are exposed to the
   *  distant star, for each coordinate direction. A negative flag value means
   *  the lower end of the box in that direction is exposed, a positive value
   *  means the upper end is exposed, and a zero value means no side is exposed
   *  in that direction. */
  char _exposed_faces[3];

  /*! @brief Box containing the simulation grid (in m). */
  Box _box;

  /*! @brief Bottom anchor of the box (in m). */
  CoordinateVector<> &_bottom_anchor;

  /*! @brief Top anchor of the box (in m). */
  CoordinateVector<> _top_anchor;

public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the distant star.
   * @param simulation_box Box containing the density grid.
   * @param log Log to write logging info to.
   */
  DistantStarContinuousPhotonSource(CoordinateVector<> position,
                                    Box simulation_box, Log *log = nullptr)
      : _position(position), _exposed_faces{0}, _box(simulation_box),
        _bottom_anchor(_box.get_anchor()), _top_anchor(_box.get_top_anchor()) {
    unsigned int num_exposed = 0;
    for (unsigned int i = 0; i < 3; ++i) {
      if (_position[i] < _bottom_anchor[i]) {
        _exposed_faces[i] = -1;
      } else if (_position[i] > _top_anchor[i]) {
        _exposed_faces[i] = 1;
      }
      num_exposed += (_exposed_faces[i] != 0);
    }

    // make sure we are really dealing with a distant star
    if (num_exposed == 0) {
      cmac_error("External stellar source lies inside the simulation box. This "
                 "will not work!");
    }

    if (log) {
      log->write_status(
          "Constructed a DistantStarContinuousPhotonSource at position [",
          _position.x(), " m, ", _position.y(), " m, ", _position.z(),
          "m], illuminating a box with anchor [", _box.get_anchor().x(), " m, ",
          _box.get_anchor().y(), " m, ", _box.get_anchor().z(),
          "m] and sides [", _box.get_sides().x(), " m, ", _box.get_sides().y(),
          " m, ", _box.get_sides().z(), "m].");
      if (num_exposed == 1) {
        log->write_status("1 side of the box is exposed to the star.");
      } else {
        log->write_status(num_exposed,
                          " sides of the box are exposed to the star.");
      }
    }
  }

  /**
   * @brief ParameterFile constructor
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  DistantStarContinuousPhotonSource(ParameterFile &params, Log *log = nullptr)
      : DistantStarContinuousPhotonSource(
            params.get_physical_vector< QUANTITY_LENGTH >(
                "continuousphotonsource:position"),
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides")),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DistantStarContinuousPhotonSource() {}

  /**
   * @brief Get the number of sides of the simulation box that is exposed to
   * radiation originating from the distant star.
   *
   * This is at least 1 and at most 3. Function should only be used for
   * debugging purposes, as there is no access to this function through the
   * general ContinuousPhotonSource interface.
   *
   * @return Number of sides of the simulation box exposed to the star.
   */
  unsigned int get_num_sides_exposed() const {
    unsigned int num_exposed = 0;
    for (unsigned int i = 0; i < 3; ++i) {
      num_exposed += (_exposed_faces[i] != 0);
    }
    return num_exposed;
  }

  /**
   * @brief Check if a photon with the given travel direction will enter the
   * simulation box at some point.
   *
   * @param direction Travel direction of the photon.
   * @param face_position Variable to store the entrance point of the photon in
   * (in m).
   * @return True if the photon will enter the box.
   */
  inline bool enters_box(CoordinateVector<> direction,
                         CoordinateVector<> &face_position) const {
    for (unsigned int i = 0; i < 3; ++i) {
      // the condition below excludes photons moving in a direction away from
      // the box
      // it also automatically excludes non-exposed faces and photons moving in
      // a plane parallel to the faces in that coordinate direction
      if (_exposed_faces[i] * direction[i] < 0.) {
        // find the intersection point
        double plane_value;
        if (_exposed_faces[i] < 0) {
          plane_value = _bottom_anchor[i];
        } else {
          plane_value = _top_anchor[i];
        }
        double l = (plane_value - _position[i]) / direction[i];
        face_position = _position + l * direction;
        unsigned int j1 = (i + 1) % 3;
        unsigned int j2 = (i + 2) % 3;
        if (face_position[j1] >= _bottom_anchor[j1] &&
            face_position[j1] <= _top_anchor[j1] &&
            face_position[j2] >= _bottom_anchor[j2] &&
            face_position[j2] <= _top_anchor[j2]) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * @brief Get the entrance point and direction of a random incoming photon.
   *
   * Since it is very hard to constrain the solid angle covered by the receiving
   * part of the box, we use a rejection technique to sample a random outgoing
   * direction from the star that ends up inside the box. The incoming direction
   * is already given, while the entrance point on the surface of the box is
   * calculated as part of the rejection technique.
   *
   * @param random_generator RandomGenerator used to generate random numbers.
   * @return std::pair containing the position and direction of a random
   * incoming photon.
   */
  virtual std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const {
    CoordinateVector<> direction =
        PhotonSource::get_random_direction(random_generator);
    CoordinateVector<> position;

    // check if the direction will enter the box at some point
    while (!enters_box(direction, position)) {
      // if not, resample
      direction = PhotonSource::get_random_direction(random_generator);
    }

    return std::make_pair(position, direction);
  }

  /**
   * @brief Get the total surface area through which radiation enters the
   * simulation box.
   *
   * This is the total area of all sides of the box that are exposed to the
   * star.
   *
   * There are at most 3 sides that can be exposed to the star, one in every
   * coordinate direction. To check whether the side in a direction is exposed,
   * we need to check if the star is above, below or inside the box in that
   * direction. If the star is inside, the side in that direction is shielded
   * from the star. If it is above or below, then respectively the top or bottom
   * side is exposed. Since top and bottom side have the same surface area, we
   * do not care about which side is exposed.
   *
   * @return Total surface area through which radiation enters the box (in m^2).
   */
  virtual double get_total_surface_area() const {
    double surface_area = 0.;

    if (_exposed_faces[0]) {
      surface_area += _box.get_sides().y() * _box.get_sides().z();
    }

    if (_exposed_faces[1]) {
      surface_area += _box.get_sides().x() * _box.get_sides().z();
    }

    if (_exposed_faces[2]) {
      surface_area += _box.get_sides().x() * _box.get_sides().y();
    }

    return surface_area;
  }
};

#endif // DISTANTSTARCONTINUOUSPHOTONSOURCE_HPP
