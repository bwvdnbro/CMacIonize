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
 * @file FLASHSnapshotDensityFunctionModule.cpp
 *
 * @brief Python module exposure of FLASHSnapshotDensityFunction.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "FLASHSnapshotDensityFunction.hpp"
#include "DensityValues.hpp"
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>

/**
 * @brief Get the temperature at a given position.
 *
 * @param density_function FLASHSnapshotDensityFunction on which we act.
 * @param x X-coordinate of a position (in m).
 * @param y Y-coordinate of a position (in m).
 * @param z Z-coordinate of a position (in m).
 * @return Temperature at that position (in K).
 */
static double get_temperature(FLASHSnapshotDensityFunction &density_function,
                              double x, double y, double z) {
  CoordinateVector<> position(x, y, z);
  DensityValues vals = density_function(position);
  return vals.get_temperature();
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE(libflashsnapshotdensityfunction) {
  boost::python::class_< FLASHSnapshotDensityFunction >(
      "FLASHSnapshotDensityFunction", boost::python::init< std::string >())
      .def("get_temperature", &get_temperature);
}
