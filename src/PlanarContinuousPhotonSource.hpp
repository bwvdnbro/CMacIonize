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
 * @file PlanarContinuousPhotonSource.hpp
 *
 * @brief Class used to generate an isotropic radiation field originating in
 * a plane perpendicular to one of the coordinate axes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PLANARCONTINUOUSPHOTONSOURCE_HPP
#define PLANARCONTINUOUSPHOTONSOURCE_HPP

#include "ContinuousPhotonSource.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Class used to generate an isotropic radiation field originating in
 * a plane perpendicular to one of the coordinate axes.
 */
class PlanarContinuousPhotonSource : public ContinuousPhotonSource {
private:
  /*! @brief Index of the coordinate that is kept fixed. */
  const uint_fast8_t _fixed_coordinate_index;

  /*! @brief Value for the coordinate that is kept fixed (in m). */
  const double _fixed_coordinate_value;

  /*! @brief Indices for the coordinates that are not kept fixed. */
  const uint_fast8_t _non_fixed_indices[2];

  /*! @brief Anchors for the non-fixed coordinates (in m). */
  const double _anchors[2];

  /*! @brief Sides for the non-fixed coordinates (in m). */
  const double _sides[2];

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
   * @brief Get the requested index for a non-fixed coordinate, given the fixed
   * coordinate index.
   *
   * A requested index of 0 will return the lowest non-fixed coordinate index,
   * an index of 1 will return the highest non-fixed coordinate index.
   *
   * @param fixed_index Index of the fixed coordinate (0/1/2).
   * @param index Requested index of the non-fixed coordinate (0/1).
   * @return Lowest/highest non-fixed coordinate index.
   */
  inline static uint_fast8_t get_non_fixed_index(const uint_fast8_t fixed_index,
                                                 const uint_fast8_t index) {

    const uint_fast8_t i0 = (fixed_index + 1) % 3;
    const uint_fast8_t i1 = (fixed_index + 2) % 3;
    const uint_fast8_t options[2] = {std::min(i0, i1), std::max(i0, i1)};
    return options[index];
  }

public:
  /**
   * @brief Constructor.
   *
   * @param fixed_coordinate_name Name of the coordinate axis perpendicular to
   * the plane (x/y/z).
   * @param fixed_coordinate_value Value for the fixed coordinate (in m).
   * @param anchor0 Anchor for the first non-fixed coordinate (in m).
   * @param anchor1 Anchor for the second non-fixed coordinate (in m).
   * @param sides0 Side for the first non-fixed coordinate (in m).
   * @param sides1 Side for the second non-fixed coordinate (in m).
   * @param luminosity Total ionizing luminosity (in s^-1).
   * @param log Log to write logging info to.
   */
  inline PlanarContinuousPhotonSource(const std::string fixed_coordinate_name,
                                      const double fixed_coordinate_value,
                                      const double anchor0,
                                      const double anchor1, const double sides0,
                                      const double sides1,
                                      const double luminosity,
                                      Log *log = nullptr)
      : _fixed_coordinate_index(get_coordinate_index(fixed_coordinate_name)),
        _fixed_coordinate_value(fixed_coordinate_value),
        _non_fixed_indices{get_non_fixed_index(_fixed_coordinate_index, 0),
                           get_non_fixed_index(_fixed_coordinate_index, 1)},
        _anchors{anchor0, anchor1}, _sides{sides0, sides1},
        _luminosity(luminosity) {

    if (log) {
      log->write_status(
          "Constructed PlanarContinuousPhotonSource perpendicular to the axis "
          "c[",
          _fixed_coordinate_index, "] = ", _fixed_coordinate_value,
          ", with bounds [", _anchors[0], ", ", _anchors[1], "] m - [",
          _sides[0], ", ", _sides[1], "m for indices [", _non_fixed_indices[0],
          ", ", _non_fixed_indices[1],
          "] and with a total ionizing luminosity of ", _luminosity, " s^-1.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the file:
   *  - normal axis: Name of the coordinate axis normal to the plane (x/y/z,
   *    default: z)
   *  - intercept: Position of the intercept along the normal axis (default:
   *    0. m)
   *  - anchor 0: Anchor for the first non-fixed coordinate (default: 0. m)
   *  - anchor 1: Anchor for the second non-fixed coordinate (default: 0. m)
   *  - side 0: Side for the first non-fixed coordinate (default: 1. m)
   *  - side 1: Side for the second non-fixed coordinate (default: 1. m)
   *  - luminosity: Total ionizing luminosity (default: 1.e48 s^-1)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline PlanarContinuousPhotonSource(ParameterFile &params, Log *log = nullptr)
      : PlanarContinuousPhotonSource(
            params.get_value< std::string >(
                "ContinuousPhotonSource:normal axis", "z"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:intercept", "0. m"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:anchor 0", "0. m"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:anchor 1", "0. m"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:side 0", "1. m"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:side 1", "1. m"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "ContinuousPhotonSource:luminosity", "1.e48 s^-1"),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~PlanarContinuousPhotonSource() {}

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
    position[_non_fixed_indices[0]] =
        _anchors[0] + random_generator.get_uniform_random_double() * _sides[0];
    position[_non_fixed_indices[1]] =
        _anchors[1] + random_generator.get_uniform_random_double() * _sides[1];
    position[_fixed_coordinate_index] = _fixed_coordinate_value;

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
    return _sides[0] * _sides[1];
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

#endif // PLANARCONTINUOUSPHOTONSOURCE_HPP
