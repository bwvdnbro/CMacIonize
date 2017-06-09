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

class VoronoiGeneratorDistribution;

/**
 * @brief DensityGrid implementation that uses an unstructured Voronoi grid.
 */
class VoronoiDensityGrid : public DensityGrid {
private:
  /*! @brief VoronoiGeneratorDistribution used to generate generator
   *  positions. */
  VoronoiGeneratorDistribution *_position_generator;

  /*! @brief Underlying Voronoi grid. */
  VoronoiGrid _voronoi_grid;

  /*! @brief Number of Lloyd iterations to apply to the grid after it has been
   *  constructed for the first time. */
  unsigned char _num_lloyd;

  /*! @brief Velocity of the grid generators (in m s^-1). */
  std::vector< CoordinateVector<> > _hydro_generator_velocity;

  /*! @brief Time step used in the hydro scheme (in s). */
  double _hydro_timestep;

  /*! @brief Polytropic index for the ideal gas equation of state. */
  double _hydro_gamma;

  /*! @brief Epsilon displacement factor used to guarantee a point lies inside
   *  a cell. */
  double _epsilon;

public:
  VoronoiDensityGrid(
      VoronoiGeneratorDistribution *position_generator,
      DensityFunction &density_function, Box box, unsigned char num_lloyd = 0,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, double hydro_timestep = 0.,
      double hydro_gamma = 5. / 3., Log *log = nullptr);

  VoronoiDensityGrid(ParameterFile &params, DensityFunction &density_function,
                     Log *log = nullptr);

  virtual ~VoronoiDensityGrid();

  virtual void initialize(std::pair< unsigned long, unsigned long > &block);
  virtual void evolve(double timestep);
  virtual void set_grid_velocity();

  virtual CoordinateVector<>
  get_interface_velocity(const iterator left, const iterator right,
                         const CoordinateVector<> interface_midpoint) const;

  virtual unsigned int get_number_of_cells() const;
  virtual unsigned long get_cell_index(CoordinateVector<> position) const;
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) const;
  virtual std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                                   CoordinateVector<>, double > >
  get_neighbours(unsigned long index);
  virtual double get_cell_volume(unsigned long index) const;
  virtual DensityGrid::iterator interact(Photon &photon, double optical_depth);
  virtual double get_total_emission(CoordinateVector<> origin,
                                    CoordinateVector<> direction,
                                    EmissionLine line);
  virtual DensityGrid::iterator begin();
  virtual DensityGrid::iterator end();
};

#endif // VORONOIDENSITYGRID_HPP
