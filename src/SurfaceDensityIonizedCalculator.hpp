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
 * @file SurfaceDensityIonizedCalculator.hpp
 *
 * @brief Object used to calculate the surface density of
 * neutral and H-alpha for a distributed grid
 * at runtime.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 *
 */
#ifndef SurfaceDensityIonizedCALCULATOR_HPP
#define SurfaceDensityIonizedCALCULATOR_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "HydroDensitySubGrid.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Object used to store the surface density values for a single subgrid.
 */
class SurfaceDensityIonizedCalculatorValues {
private:
  /*! @brief Surface density values (in kg m^-2). */
  std::vector< double > _surface_densities_Halpha;
  std::vector< double > _surface_densities_neutral;

  /*! @brief Number of pixels in the vertical direction. */
  const int_fast32_t _ny;

public:
  /**
   * @brief Constructor.
   *
   * @param nx Number of horizontal pixels.
   * @param ny Number of vertical pixels.
   */
  inline SurfaceDensityIonizedCalculatorValues(const int_fast32_t nx,
                                               const int_fast32_t ny)
      : _surface_densities_Halpha(nx * ny, 0.),
        _surface_densities_neutral(nx * ny, 0.), _ny(ny) {}

  /**
   * @brief Reset the pixel values to zero.
   */
  inline void reset() {
    for (uint_fast32_t i = 0; i < _surface_densities_neutral.size(); ++i) {
      _surface_densities_neutral[i] = 0.;
      _surface_densities_Halpha[i] = 0.;
    }
  }

  /**
   * @brief Access the pixel with the given index values.
   *
   * @param ix Horizontal index.
   * @param iy Vertical index.
   * @return Reference to the corresponding pixel.
   */
  inline double &pixel_neutral(const int_fast32_t ix, const int_fast32_t iy) {
    return _surface_densities_neutral[(ix * _ny + iy)];
  }
  inline double &pixel_ion(const int_fast32_t ix, const int_fast32_t iy) {
    return _surface_densities_Halpha[(ix * _ny + iy)];
  }

  /**
   * @brief Multiply all pixels with the given factor.
   *
   * @param factor Factor.
   * @return Reference to the updated object.
   */
  inline SurfaceDensityIonizedCalculatorValues &
  operator*=(const double factor) {
    for (uint_fast32_t i = 0; i < _surface_densities_neutral.size(); ++i) {
      _surface_densities_neutral[i] *= factor;
      _surface_densities_Halpha[i] *= factor;
    }
    return *this;
  }

  /**
   * @brief Add the given surface density map to this one.
   *
   * @param other Other SurfaceDensityIonizedCalculatorValues object to add.
   * @return Reference to the updated object.
   */
  inline SurfaceDensityIonizedCalculatorValues &
  operator+=(const SurfaceDensityIonizedCalculatorValues &other) {
    for (uint_fast32_t i = 0; i < _surface_densities_neutral.size(); ++i) {
      _surface_densities_neutral[i] += other._surface_densities_neutral[i];
      _surface_densities_Halpha[i] += other._surface_densities_Halpha[i];
    }
    return *this;
  }
};

/**
 * @brief Object used to calculate the surface density for a distributed grid
 * at runtime.
 */
class SurfaceDensityIonizedCalculator {
private:
  /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< int_fast32_t > _number_of_subgrids;

  /*! @brief Number of cells per coordinate direction for a single subgrid. */
  const CoordinateVector< int_fast32_t > _number_of_cells;

  /*! @brief Per subgrid surface density values. */
  std::vector< SurfaceDensityIonizedCalculatorValues * > _subgrid_values;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param number_of_cells Number of cells per coordinate direction for a
   * single subgrid.
   * @param number_of_cells
   */
  inline SurfaceDensityIonizedCalculator(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells)
      : _number_of_subgrids(number_of_subgrids),
        _number_of_cells(number_of_cells),
        _subgrid_values(number_of_subgrids.x() * number_of_subgrids.y() *
                            number_of_subgrids.z(),
                        nullptr) {
    for (uint_fast32_t i = 0; i < _subgrid_values.size(); ++i) {
      _subgrid_values[i] = new SurfaceDensityIonizedCalculatorValues(
          number_of_cells.x(), number_of_cells.y());
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~SurfaceDensityIonizedCalculator() {
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
        double neutral = 0.;
        double ion = 0.;
        for (int_fast32_t iz = 0; iz < _number_of_cells.z(); ++iz) {
          const double N_fraction =
              HydroDensitySubGrid::hydroiterator(cell_index_xy + iz, subgrid)
                  .get_ionization_variables()
                  .get_ionic_fraction(ION_H_n);
          const double Density =
              HydroDensitySubGrid::hydroiterator(cell_index_xy + iz, subgrid)
                  .get_hydro_variables()
                  .get_primitives_density();

          neutral += Density * N_fraction;

          const double H_alpha = Density * (1 - N_fraction);

          ion += H_alpha * H_alpha;
        }
        _subgrid_values[index]->pixel_neutral(ix, iy) = neutral;
        _subgrid_values[index]->pixel_ion(ix, iy) = ion;
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
        SurfaceDensityIonizedCalculatorValues values(_number_of_cells.x(),
                                                     _number_of_cells.y());
        for (int_fast32_t iz = 0; iz < _number_of_subgrids.z(); ++iz) {
          const int_fast32_t subgrid_index = subgrid_index_xy + iz;
          values += *_subgrid_values[subgrid_index];
        }
        values *= box.get_sides().z();
        for (int_fast32_t cix = 0; cix < _number_of_cells.x(); ++cix) {
          for (int_fast32_t ciy = 0; ciy < _number_of_cells.y(); ++ciy) {
            file << values.pixel_neutral(cix, ciy) << "\t"
                 << values.pixel_ion(cix, ciy) << "\n";
          }
        }
      }
    }
  }
};

#endif // SurfaceDensityIonizedCALCULATOR_HPP
