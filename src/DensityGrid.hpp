/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityGrid.hpp
 *
 * @brief Density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRID_HPP
#define DENSITYGRID_HPP

#include "Box.hpp"

class CrossSections;
class DensityFunction;
class DensityValues;
class ParameterFile;
class Photon;
class RecombinationRates;

/**
 * @brief Density grid.
 *
 * Contains the actual cells with densities and neutral fractions, and the
 * routines used to calculate the optical depth along a photon path.
 */
class DensityGrid {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector<> _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells per dimension. Note that by choosing an unsigned
   *  char type, we automatically limit grid sizes to 256^3. */
  CoordinateVector< unsigned char > _ncell;

  /*! @brief Helium abundance. */
  double _helium_abundance;

  /*! @brief Density grid. */
  DensityValues ***_density;

  /*! @brief Cross sections for photoionization. */
  CrossSections &_cross_sections;

  /*! @brief Recombination rates used in ionization balance calculation. */
  RecombinationRates &_recombination_rates;

public:
  DensityGrid(Box box, CoordinateVector< unsigned char > ncell,
              double helium_abundance, double initial_temperature,
              DensityFunction &density_function, CrossSections &cross_sections,
              RecombinationRates &recombination_rates);

  DensityGrid(ParameterFile &parameters, Box box,
              CoordinateVector< unsigned char > ncell,
              DensityFunction &density_function, CrossSections &cross_sections,
              RecombinationRates &recombination_rates);

  ~DensityGrid();

  double get_total_mass();

  CoordinateVector< int > get_cell_indices(CoordinateVector<> position);
  Box get_cell(CoordinateVector< int > index);
  bool is_inside(CoordinateVector< int > index);
  CoordinateVector<> get_wall_intersection(CoordinateVector<> &photon_origin,
                                           CoordinateVector<> &photon_direction,
                                           Box &cell,
                                           CoordinateVector< char > &next_index,
                                           double &ds);

  bool interact(Photon &photon, double optical_depth);

  void find_H0(double ch1, double ch2, double che, double AHe, double T,
               double &h0, double &he0);

  void calculate_ionization_state(unsigned int nphoton);
};

#endif // DENSITYGRID_HPP
