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
 * @file GadgetSnapshotDensityFunction.cpp
 *
 * @brief DensityFunction that reads a density field from a Gadget HDF5 snapshot
 * file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"
using namespace std;

/**
 * @brief Cubic spline kernel used in Gadget2.
 *
 * @param u Distance in units of the smoothing length.
 * @param h Smoothing length.
 * @return Value of the cubic spline kernel.
 */
double GadgetSnapshotDensityFunction::cubic_spline_kernel(double u, double h) {
  const double KC1 = 2.546479089470;
  const double KC2 = 15.278874536822;
  const double KC5 = 5.092958178941;
  if (u < 1.) {
    if (u < 0.5) {
      return (KC1 + KC2 * (u - 1.) * u * u) / (h * h * h);
    } else {
      return KC5 * (1. - u) * (1. - u) * (1. - u) / (h * h * h);
    }
  } else {
    // the cubic spline kernel has compact support
    return 0.;
  }
}

/**
 * @brief Constructor.
 *
 * @param name Name of the snapshot file to read.
 */
GadgetSnapshotDensityFunction::GadgetSnapshotDensityFunction(std::string name) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);
  // open the group containing the SPH particle data
  HDF5Tools::HDF5Group gasparticles = HDF5Tools::open_group(file, "/PartType0");
  // read the positions, masses and smoothing lengths
  _positions =
      HDF5Tools::read_dataset<CoordinateVector>(gasparticles, "Coordinates");
  _masses = HDF5Tools::read_dataset<double>(gasparticles, "Masses");
  _smoothing_lengths =
      HDF5Tools::read_dataset<double>(gasparticles, "SmoothingLength");
  // close the group
  HDF5Tools::close_group(gasparticles);
  // close the file
  HDF5Tools::close_file(file);
}

/***
 * @brief Function that returns the density for the given coordinate.
 *
 * @param position CoordinateVector specifying a coordinate position.
 * @return Density at the given coordinate.
 */
double GadgetSnapshotDensityFunction::operator()(CoordinateVector position) {
  return 0.;
}
