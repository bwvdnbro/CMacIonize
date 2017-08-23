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

class SimulationBox;
class VoronoiGeneratorDistribution;
class VoronoiGrid;

/**
 * @brief DensityGrid implementation that uses an unstructured Voronoi grid.
 */
class VoronoiDensityGrid : public DensityGrid {
private:
  /*! @brief VoronoiGeneratorDistribution used to generate generator
   *  positions. */
  VoronoiGeneratorDistribution *_position_generator;

  /*! @brief Underlying Voronoi grid. */
  VoronoiGrid *_voronoi_grid;

  /*! @brief Periodicity flags for the simulation box. */
  CoordinateVector< bool > _periodicity_flags;

  /*! @brief Number of Lloyd iterations to apply to the grid after it has been
   *  constructed for the first time. */
  unsigned char _num_lloyd;

  /*! @brief Generator positions (in m). */
  std::vector< CoordinateVector<> > _generator_positions;

  /*! @brief Velocity of the grid generators (in m s^-1). */
  std::vector< CoordinateVector<> > _hydro_generator_velocity;

  /*! @brief Epsilon displacement factor used to guarantee a point lies inside
   *  a cell. */
  double _epsilon;

  /*! @brief Type of Voronoi grid to use. */
  std::string _voronoi_grid_type;

public:
  VoronoiDensityGrid(
      VoronoiGeneratorDistribution *position_generator,
      const Box<> &simulation_box, std::string grid_type = "Old",
      unsigned char num_lloyd = 0,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr);

  VoronoiDensityGrid(const SimulationBox &simulation_box, ParameterFile &params,
                     bool hydro = false, Log *log = nullptr);

  virtual ~VoronoiDensityGrid();

  virtual void initialize(std::pair< unsigned long, unsigned long > &block,
                          DensityFunction &density_function);
  virtual void evolve(double timestep);
  virtual void set_grid_velocity(double gamma);

  virtual CoordinateVector<>
  get_interface_velocity(const iterator left, const iterator right,
                         const CoordinateVector<> interface_midpoint) const;

  virtual unsigned int get_number_of_cells() const;
  virtual unsigned long get_cell_index(CoordinateVector<> position) const;
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) const;
  virtual std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                                   CoordinateVector<>, double > >
  get_neighbours(unsigned long index);
  virtual std::vector< Face > get_faces(unsigned long index) const;
  virtual double get_cell_volume(unsigned long index) const;
  virtual double integrate_optical_depth(const Photon &photon);
  virtual DensityGrid::iterator interact(Photon &photon, double optical_depth);
  virtual double get_total_emission(CoordinateVector<> origin,
                                    CoordinateVector<> direction,
                                    EmissionLine line);
  virtual DensityGrid::iterator begin();
  virtual DensityGrid::iterator end();
};

#endif // VORONOIDENSITYGRID_HPP
