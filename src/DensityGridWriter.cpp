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
 * @file DensityGridWriter.cpp
 *
 * @brief DensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DensityGridWriter.hpp"
#include "DensityGrid.hpp"
#include "HDF5Tools.hpp"

/**
 * @brief Constructor.
 *
 * @param name Name of the file to write.
 * @param grid DensityGrid containing the data to write.
 */
DensityGridWriter::DensityGridWriter(std::string name, DensityGrid &grid)
    : _name(name), _grid(grid) {}

/**
 * @brief Write the file.
 */
void DensityGridWriter::write() {
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(_name, HDF5Tools::HDF5FILEMODE_WRITE);

  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  CoordinateVector<> box(1.);
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", box);
  int dimension = 3;
  HDF5Tools::write_attribute< int >(group, "Dimension", dimension);
}
