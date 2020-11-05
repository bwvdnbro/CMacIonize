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
 * @file SurfaceDensityCalculator.hpp
 *
 * @brief Object used to calculate the surface density for a distributed grid
 * at runtime.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SURFACEDENSITYCALCULATOR_HPP
#define SURFACEDENSITYCALCULATOR_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "HydroDensitySubGrid.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Object used to store the surface density values for a single subgrid.
 */
class SurfaceDensityCalculatorValues {
private:
  /*! @brief Surface density values (in kg m^-2). */
  std::vector< double > _surface_densities;

  /*! @brief Number of pixels in the vertical direction. */
  const int_fast32_t _ny;

public:
  /**
   * @brief Constructor.
   *
   * @param nx Number of horizontal pixels.
   * @param ny Number of vertical pixels.
   */
  inline SurfaceDensityCalculatorValues(const int_fast32_t nx,
                                        const int_fast32_t ny)
      : _surface_densities(nx * ny, 0.), _ny(ny) {}

  /**
   * @brief Reset the pixel values to zero.
   */
  inline void reset() {
    for (uint_fast32_t i = 0; i < _surface_densities.size(); ++i) {
      _surface_densities[i] = 0.;
    }
  }

  /**
   * @brief Access the pixel with the given index values.
   *
   * @param ix Horizontal index.
   * @param iy Vertical index.
   * @return Reference to the corresponding pixel.
   */
  inline double &pixel(const int_fast32_t ix, const int_fast32_t iy) {
    return _surface_densities[ix * _ny + iy];
  }

  /**
   * @brief Multiply all pixels with the given factor.
   *
   * @param factor Factor.
   * @return Reference to the updated object.
   */
  inline SurfaceDensityCalculatorValues &operator*=(const double factor) {
    for (uint_fast32_t i = 0; i < _surface_densities.size(); ++i) {
      _surface_densities[i] *= factor;
    }
    return *this;
  }

  /**
   * @brief Add the given surface density map to this one.
   *
   * @param other Other SurfaceDensityCalculatorValues object to add.
   * @return Reference to the updated object.
   */
  inline SurfaceDensityCalculatorValues &
  operator+=(const SurfaceDensityCalculatorValues &other) {
    for (uint_fast32_t i = 0; i < _surface_densities.size(); ++i) {
      _surface_densities[i] += other._surface_densities[i];
    }
    return *this;
  }
};

/**
 * @brief Object used to calculate the surface density for a distributed grid
 * at runtime.
 */
class SurfaceDensityCalculator {
private:
  /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< int_fast32_t > _number_of_subgrids;

  /*! @brief Number of cells per coordinate direction for a single subgrid. */
  const CoordinateVector< int_fast32_t > _number_of_cells;

  /*! @brief Per subgrid surface density values. */
  std::vector< SurfaceDensityCalculatorValues * > _subgrid_values;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param number_of_cells Number of cells per coordinate direction for a
   * single subgrid.
   */
  inline SurfaceDensityCalculator(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells)
      : _number_of_subgrids(number_of_subgrids),
        _number_of_cells(number_of_cells),
        _subgrid_values(number_of_subgrids.x() * number_of_subgrids.y() *
                            number_of_subgrids.z(),
                        nullptr) {
    for (uint_fast32_t i = 0; i < _subgrid_values.size(); ++i) {
      _subgrid_values[i] = new SurfaceDensityCalculatorValues(
          number_of_cells.x(), number_of_cells.y());
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~SurfaceDensityCalculator() {
    for (uint_fast32_t i = 0; i < _subgrid_values.size(); ++i) {
      delete _subgrid_values[i];
    }
  }

  /**
   * @brief Calculate the surface densities for the given subgrid.
   *
   * @param index Index of the subgrid.
   * @param subgrid Subgrid.
   */
  inline void calculate_surface_density(const uint_fast32_t index,
                                        HydroDensitySubGrid &subgrid) {

    _subgrid_values[index]->reset();
    for (int_fast32_t ix = 0; ix < _number_of_cells.x(); ++ix) {
      const int_fast32_t cell_index_x =
          ix * _number_of_cells.y() * _number_of_cells.z();
      for (int_fast32_t iy = 0; iy < _number_of_cells.y(); ++iy) {
        const int_fast32_t cell_index_xy =
            cell_index_x + iy * _number_of_cells.z();
        double value = 0.;
        for (int_fast32_t iz = 0; iz < _number_of_cells.z(); ++iz) {
          value +=
              HydroDensitySubGrid::hydroiterator(cell_index_xy + iz, subgrid)
                  .get_hydro_variables()
                  .get_primitives_density();
        }
        _subgrid_values[index]->pixel(ix, iy) = value;
      }
    }
  }

  /**
   * @brief Output the surface densities.
   *
   * @param filename Name of the file to write.
   * @param box Size of the simulation box (in m).
   */
  inline void output(const std::string filename, const Box<> box) const {

    std::ofstream file(filename);
    file << _number_of_subgrids.x() << "\t" << _number_of_subgrids.y() << "\n";
    file << _number_of_cells.x() << "\t" << _number_of_cells.y() << "\n";
    file << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\n";
    file << box.get_sides().x() << "\t" << box.get_sides().y() << "\n";
    for (int_fast32_t ix = 0; ix < _number_of_subgrids.x(); ++ix) {
      const int_fast32_t subgrid_index_x =
          ix * _number_of_subgrids.y() * _number_of_subgrids.z();
      for (int_fast32_t iy = 0; iy < _number_of_subgrids.y(); ++iy) {
        const int_fast32_t subgrid_index_xy =
            subgrid_index_x + iy * _number_of_subgrids.z();
        SurfaceDensityCalculatorValues values(_number_of_cells.x(),
                                              _number_of_cells.y());
        for (int_fast32_t iz = 0; iz < _number_of_subgrids.z(); ++iz) {
          const int_fast32_t subgrid_index = subgrid_index_xy + iz;
          values += *_subgrid_values[subgrid_index];
        }
        values *= box.get_sides().z();
        for (int_fast32_t cix = 0; cix < _number_of_cells.x(); ++cix) {
          for (int_fast32_t ciy = 0; ciy < _number_of_cells.y(); ++ciy) {
            file << values.pixel(cix, ciy) << "\n";
          }
        }
      }
    }
  }
};

#endif // SURFACEDENSITYCALCULATOR_HPP
