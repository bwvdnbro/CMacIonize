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
 * @file testParallelCartesianDensityGrid.cpp
 *
 * @brief Unit test for the ParallelCartesianDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "BlockSyntaxDensityFunction.hpp"
#include "Box.hpp"
#include "HDF5Tools.hpp"
#include "ParallelCartesianDensityGrid.hpp"

/**
 * @brief Unit test for the ParallelCartesianDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // this is shorthand for box(CoordinateVector<>(0.), CoordinateVector<>(1.))
  // and only works because we have defined the single value constructor for
  // CoordinateVector<>
  Box<> box(0., 1.);
  CoordinateVector< int > numcell(32, 32, 32);
  unsigned int numdomain = 64;
  std::pair< int, int > domain = std::make_pair(0, 64);

  ParallelCartesianDensityGrid grid(box, numcell, numdomain, domain);

  BlockSyntaxDensityFunction density_function("blocksyntaxtest.yml");

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    (*it).initialize(density_function);
  }

  unsigned int totnumcell = grid.get_number_of_cells();

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file("paralleltest.hdf5", HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int dimension = 3;
  HDF5Tools::write_attribute< int >(group, "Dimension", dimension);
  std::vector< unsigned int > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int numfiles = 1;
  HDF5Tools::write_attribute< int >(group, "NumFilesPerSnapshot", numfiles);
  std::vector< unsigned int > numpart(6, 0);
  numpart[0] = totnumcell;
  std::vector< unsigned int > numpart_high(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total_HighWord", numpart_high);
  double time = 0.;
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_value = 1;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_value);
  HDF5Tools::close_group(group);

  group = HDF5Tools::create_group(file, "/PartType0");
  HDF5Tools::create_dataset< CoordinateVector<> >(group, "Coordinates",
                                                  totnumcell);
  HDF5Tools::create_dataset< double >(group, "NumberDensity", totnumcell);
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    ParallelCartesianDensitySubGrid &subgrid =
        reinterpret_cast< ParallelCartesianDensitySubGrid & >(*it);
    std::vector< CoordinateVector<> > positions = subgrid.get_positions();
    HDF5Tools::append_dataset(group, "Coordinates", it.offset(), positions);
    HDF5Tools::append_dataset(group, "NumberDensity", it.offset(),
                              subgrid.get_number_density_handle());
  }
  HDF5Tools::close_group(group);

  HDF5Tools::close_file(file);

  return 0;
}
