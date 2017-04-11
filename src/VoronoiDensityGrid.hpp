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
 * @file VoronoiDensityGrid.hpp
 *
 * @brief DensityGrid implementation that uses an unstructured Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIDENSITYGRID_HPP
#define VORONOIDENSITYGRID_HPP

#include "DensityGrid.hpp"
#include "VoronoiGrid.hpp"

class VoronoiDensityGridPositions;

/**
 * @brief DensityGrid implementation that uses an unstructured Voronoi grid.
 */
class VoronoiDensityGrid : public DensityGrid {
private:
  /*! @brief VoronoiDensityGridPositions object used to generate generator
   *  positions. */
  VoronoiDensityGridPositions &_position_generator;

  /*! @brief Underlying Voronoi grid. */
  VoronoiGrid _voronoi_grid;

public:
  VoronoiDensityGrid(
      VoronoiDensityGridPositions &position_generator,
      DensityFunction &density_function, Box box,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~VoronoiDensityGrid() {}

  virtual void initialize(std::pair< unsigned long, unsigned long > &block);

  virtual unsigned int get_number_of_cells() const;
  virtual unsigned long get_cell_index(CoordinateVector<> position) const;
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) const;
  virtual std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                                   CoordinateVector<>, double > >
  get_neighbours(unsigned long index);
  virtual double get_cell_volume(unsigned long index) const;
  virtual DensityGrid::iterator interact(Photon &photon, double optical_depth);
  virtual DensityGrid::iterator begin();
  virtual DensityGrid::iterator end();
};

#endif // VORONOIDENSITYGRID_HPP
