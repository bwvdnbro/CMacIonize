/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SPHNGVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution implementation that reads Voronoi
 * generator sites from an SPHNG binary snapshot dump.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPHNGVORONOIGENERATORDISTRIBUTION_HPP
#define SPHNGVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <vector>

/**
 * @brief VoronoiGeneratorDistribution implementation that reads Voronoi
 * generator sites from an SPHNG binary snapshot dump.
 */
class SPHNGVoronoiGeneratorDistribution : public VoronoiGeneratorDistribution {
private:
  /**
   * @brief Functor used to filter out positions that are not inside the
   * simulation box.
   */
  class BoxFilter {
  private:
    /*! @brief Underlying box (in m). */
    const Box<> &_box;

  public:
    /**
     * @brief Constructor.
     *
     * @param box Underlying box (in m).
     */
    inline BoxFilter(const Box<> &box) : _box(box) {}

    /**
     * @brief Return false if the given position is inside the box.
     *
     * @param p Position (in m).
     * @return True if the position is outside the box.
     */
    inline bool operator()(const CoordinateVector<> &p) const {
      return !_box.inside(p);
    }
  };

  /*! @brief Box containing the generators (in m). */
  const Box<> _box;

  /*! @brief Number of positions already generated. */
  generatornumber_t _current_number;

  /*! @brief Generator positions. */
  std::vector< CoordinateVector<> > _generator_positions;

public:
  /**
   * @brief Comparison function for two 3D positions.
   *
   * We will use the following rules to decide if a position
   * @$fa=(a_x,a_y,a_z)@f$ is smaller than a position @$fb=(b_x,b_y,b_z)@f$:
   *  - the position is smaller if @f$a_x<b_x@f$,
   *  - if @f$a_x=b_x@f$, then it is smaller if @f$a_y<b_y@f$,
   *  - if @f$a_y=b_y@f$, then it is smaller if @f$a_z<b_z@f$,
   *  - if @f$a_z=b_z@f$, we assume it is larger.
   *
   * @param a Position @f$a@f$.
   * @param b Position @f$b@f$.
   * @return True if @f$a<b@f$, according to the rules above.
   */
  inline static bool position_smaller_than(const CoordinateVector<> &a,
                                           const CoordinateVector<> &b) {
    if (a.x() < b.x()) {
      return true;
    } else {
      if (a.x() == b.x()) {
        if (a.y() < b.y()) {
          return true;
        } else {
          if (a.y() == b.y()) {
            if (a.z() < b.z()) {
              return true;
            } else {
              return false;
            }
          } else {
            return false;
          }
        }
      } else {
        return false;
      }
    }
  }

  SPHNGVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                    std::string filename, Log *log = nullptr);

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - filename: Name of the file (required)
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SPHNGVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                    ParameterFile &params, Log *log = nullptr)
      : SPHNGVoronoiGeneratorDistribution(
            simulation_box,
            params.get_filename(
                "DensityGrid:VoronoiGeneratorDistribution:filename"),
            log) {}

  /**
   * @brief Get the number of positions that this distribution generates.
   *
   * @return Number of generated positions.
   */
  virtual generatornumber_t get_number_of_positions() const {
    return _generator_positions.size();
  }

  /**
   * @brief Get SPH generator position.
   *
   * @return SPH generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    cmac_assert(_current_number < _generator_positions.size());
    CoordinateVector<> pos = _generator_positions[_current_number];
    ++_current_number;
    return pos;
  }
};

#endif // SPHNGVORONOIGENERATORDISTRIBUTION_HPP
