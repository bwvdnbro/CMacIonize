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
 * @file ExtendedDiscContinuousPhotonSource.hpp
 *
 * @brief Class used to generate an isotropic radiation field originating in
 * an extended Gaussian disc perpendicular to one of the coordinate axes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EXTENDEDDISCCONTINUOUSPHOTONSOURCE_HPP
#define EXTENDEDDISCCONTINUOUSPHOTONSOURCE_HPP

#include "Box.hpp"
#include "ContinuousPhotonSource.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Class used to generate an isotropic radiation field originating in
 * an extended Gaussian disc perpendicular to one of the coordinate axes.
 */
class ExtendedDiscContinuousPhotonSource : public ContinuousPhotonSource {
private:
  /*! @brief Dimensions of the simulation box (in m). */
  const Box<> _simulation_box;

  /*! @brief Index of the coordinate axis perpendicular to the disc. */
  const uint_fast8_t _disc_coordinate_index;

  /*! @brief Origin of the disc (in m). */
  const double _disc_coordinate_origin;

  /*! @brief Scale height of the disc (in m). */
  const double _disc_coordinate_scale_height;

  /*! @brief Indices for the coordinates in the plane of the disc. */
  const uint_fast8_t _disc_plane_indices[2];

  /*! @brief Total ionizing luminosity (in s^-1). */
  const double _luminosity;

  /**
   * @brief Get the index corresponding to the given coordinate axis name.
   *
   * @param name Name of a coordinate axis (x/y/z).
   * @return Corresponding index (0/1/2).
   */
  inline static uint_fast8_t get_coordinate_index(const std::string name) {

    if (name == "x") {
      return 0;
    } else if (name == "y") {
      return 1;
    } else if (name == "z") {
      return 2;
    } else {
      cmac_error("Unknown coordinate axis name: \"%s\"!", name.c_str());
    }
  }

  /**
   * @brief Get the requested index for a disc plane coordinate, given the index
   * of the coordinate axis perpendicular to the plane.
   *
   * A requested index of 0 will return the lowest disc coordinate index,
   * an index of 1 will return the highest disc coordinate index.
   *
   * @param disc_index Index of the coordinate axis perpendicular to the disc
   * (0/1/2).
   * @param index Requested index of the disc coordinate (0/1).
   * @return Lowest/highest disc coordinate index.
   */
  inline static uint_fast8_t get_disc_plane_index(const uint_fast8_t disc_index,
                                                  const uint_fast8_t index) {

    const uint_fast8_t i0 = (disc_index + 1) % 3;
    const uint_fast8_t i1 = (disc_index + 2) % 3;
    const uint_fast8_t options[2] = {std::min(i0, i1), std::max(i0, i1)};
    return options[index];
  }

public:
  /**
   * @brief Constructor.
   *
   * @param simulation_box Dimensions of the simulation box (in m).
   * @param disc_coordinate_name Name of the coordinate axis perpendicular to
   * the plane of the disc (x/y/z).
   * @param disc_coordinate_origin Origin of the disc (in m).
   * @param disc_coordinate_scale_height Scale height of the disc (in m).
   * @param luminosity Total ionizing luminosity (in s^-1).
   * @param log Log to write logging info to.
   */
  inline ExtendedDiscContinuousPhotonSource(
      const Box<> &simulation_box, const std::string disc_coordinate_name,
      const double disc_coordinate_origin,
      const double disc_coordinate_scale_height, const double luminosity,
      Log *log = nullptr)
      : _simulation_box(simulation_box),
        _disc_coordinate_index(get_coordinate_index(disc_coordinate_name)),
        _disc_coordinate_origin(disc_coordinate_origin),
        _disc_coordinate_scale_height(disc_coordinate_scale_height),
        _disc_plane_indices{get_disc_plane_index(_disc_coordinate_index, 0),
                            get_disc_plane_index(_disc_coordinate_index, 1)},
        _luminosity(luminosity) {

    if (log) {
      log->write_status(
          "Constructed PlanarContinuousPhotonSource perpendicular to the axis "
          "c[",
          _disc_coordinate_index, "], with bounds [",
          _simulation_box.get_anchor()[_disc_plane_indices[0]], ", ",
          _simulation_box.get_anchor()[_disc_plane_indices[1]], "] m - [",
          _simulation_box.get_sides()[_disc_plane_indices[0]], ", ",
          _simulation_box.get_sides()[_disc_plane_indices[1]],
          "m for indices [", _disc_plane_indices[0], ", ",
          _disc_plane_indices[1], "], origin ", _disc_coordinate_origin,
          " m, scale height ", _disc_coordinate_scale_height,
          " m and with a total ionizing luminosity of ", _luminosity, " s^-1.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the file:
   *  - normal axis: Name of the coordinate axis normal to the plane of the disc
   *    (x/y/z, default: z)
   *  - intercept: Position of the intercept along the normal axis, which
   *    defines the origin of the disc (default: 0. m)
   *  - scale height: Scale height of the disc (default: 200. pc)
   *  - luminosity: Total ionizing luminosity (default: 1.e48 s^-1)
   *
   * @param simulation_box Dimensions of the simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline ExtendedDiscContinuousPhotonSource(const Box<> simulation_box,
                                            ParameterFile &params,
                                            Log *log = nullptr)
      : ExtendedDiscContinuousPhotonSource(
            simulation_box,
            params.get_value< std::string >(
                "ContinuousPhotonSource:normal axis", "z"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:intercept", "0. m"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:scale height", "200. pc"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "ContinuousPhotonSource:luminosity", "1.e48 s^-1"),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~ExtendedDiscContinuousPhotonSource() {}

  /**
   * @brief Get the starting position and direction of a random planar photon.
   *
   * @param random_generator RandomGenerator to use.
   * @return std::pair of CoordinateVector instances, specifying a starting
   * position and direction for a photon.
   */
  virtual std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const {

    CoordinateVector<> position;
    // disc plane
    position[_disc_plane_indices[0]] =
        _simulation_box.get_anchor()[_disc_plane_indices[0]] +
        random_generator.get_uniform_random_double() *
            _simulation_box.get_sides()[_disc_plane_indices[0]];
    position[_disc_plane_indices[1]] =
        _simulation_box.get_anchor()[_disc_plane_indices[1]] +
        random_generator.get_uniform_random_double() *
            _simulation_box.get_sides()[_disc_plane_indices[1]];

    // disc extent
    position[_disc_coordinate_index] =
        _disc_coordinate_scale_height *
            std::sqrt(-2. *
                      std::log(random_generator.get_uniform_random_double())) *
            std::cos(2. * M_PI * random_generator.get_uniform_random_double()) +
        _disc_coordinate_origin;
    while (position[_disc_coordinate_index] <
               _simulation_box.get_anchor()[_disc_coordinate_index] ||
           position[_disc_coordinate_index] >
               _simulation_box.get_anchor()[_disc_coordinate_index] +
                   _simulation_box.get_sides()[_disc_coordinate_index]) {
      position[_disc_coordinate_index] =
          _disc_coordinate_scale_height *
              std::sqrt(
                  -2. *
                  std::log(random_generator.get_uniform_random_double())) *
              std::cos(2. * M_PI *
                       random_generator.get_uniform_random_double()) +
          _disc_coordinate_origin;
    }

    // random direction
    const double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(0., 1. - cost * cost));
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);
    const CoordinateVector<> direction(sint * cosp, sint * sinp, cost);

    return std::make_pair(position, direction);
  }

  /**
   * @brief Get the total surface area through which the radiation enters the
   * simulation box.
   *
   * @return Total surface area (in m^2).
   */
  virtual double get_total_surface_area() const {
    cmac_error("This function should not be used!");
    return 0.;
  }

  /**
   * @brief Does this ContinuousPhotonSource have a total luminosity value?
   *
   * @return True.
   */
  virtual bool has_total_luminosity() const { return true; }

  /**
   * @brief Get the total luminosity for the source.
   *
   * @return Total ionizing luminosity of the source (in s^-1).
   */
  virtual double get_total_luminosity() const { return _luminosity; }
};

#endif // EXTENDEDDISCCONTINUOUSPHOTONSOURCE_HPP
