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
 * @file NewVoronoiCell.hpp
 *
 * @brief Voronoi cell constructed using the new algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOICELL_HPP
#define NEWVORONOICELL_HPP

#include "VoronoiFace.hpp"

#include <vector>

/**
 * @brief Voronoi cell constructed using the new algorithm.
 */
class NewVoronoiCell {
private:
  /*! @brief Volume of the cell (in m^3). */
  double _volume;

  /*! @brief Centroid of the cell (in m). */
  CoordinateVector<> _centroid;

  /*! @brief Faces of the cell. */
  std::vector< VoronoiFace > _faces;

public:
  /**
   * @brief Empty constructor.
   */
  inline NewVoronoiCell() : _volume(0.) {}

  /**
   * @brief Constructor.
   *
   * @param volume Volume of the cell (in m^3).
   * @param centroid Centroid of the cell (in m).
   * @param faces Faces of the cell.
   */
  inline NewVoronoiCell(double volume, const CoordinateVector<> &centroid,
                        const std::vector< VoronoiFace > &faces)
      : _volume(volume), _centroid(centroid), _faces(faces) {}

  /**
   * @brief Get the volume of the cell.
   *
   * @return Volume of the cell (in m^3).
   */
  inline double get_volume() const { return _volume; }

  /**
   * @brief Get the centroid of the cell.
   *
   * @return Centroid of the cell (in m).
   */
  inline const CoordinateVector<> &get_centroid() const { return _centroid; }

  /**
   * @brief Get the faces of the cell.
   *
   * @return Faces of the cell.
   */
  inline const std::vector< VoronoiFace > &get_faces() const { return _faces; }
};

#endif // NEWVORONOICELL_HPP
