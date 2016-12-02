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
 * @file DensityGridModule.cpp
 *
 * @brief Python module exposure of DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CMacIonizeSnapshotDensityFunction.hpp"
#include "CartesianDensityGrid.hpp"
#include "DensityGridFactory.hpp"
#include "HDF5Tools.hpp"
#include "ParameterFile.hpp"
#include <boost/noncopyable.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/module.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

/**
 * @brief Python constructor for a DensityGrid.
 *
 * Python does not know about different types of DensityGrid. We just provide a
 * single interface, with all type specific stuff being handled by C++. A
 * Python DensityGrid is constructed from a snapshot file, which contains all
 * the relevant info we need to figure out which type we are actually using.
 *
 * @param filename Name of the snapshot file containing the grid.
 * @return boost::shared_ptr to the DensityGrid.
 */
static boost::shared_ptr< DensityGrid >
initDensityGrid(const std::string &filename) {
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // read parameters
  HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/Parameters");
  std::vector< std::string > parameternames =
      HDF5Tools::get_attribute_names(group);
  ParameterFile parameters;
  for (auto it = parameternames.begin(); it != parameternames.end(); ++it) {
    std::string attname = *it;
    std::string attvalue =
        HDF5Tools::read_attribute< std::string >(group, attname);
    parameters.add_value(attname, attvalue);
  }
  HDF5Tools::close_group(group);

  HDF5Tools::close_file(file);

  CMacIonizeSnapshotDensityFunction density_function(filename);

  return boost::shared_ptr< DensityGrid >(
      DensityGridFactory::generate(parameters, density_function));
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE(libdensitygrid) {
  // we need to use no_init to tell Boost that we provide a custom constructor.
  // we need to use noncopyable to tell Boost that we want to construct an
  // object of an abstract type.
  boost::python::class_< DensityGrid, boost::shared_ptr< DensityGrid >,
                         boost::noncopyable >("DensityGrid",
                                              boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&initDensityGrid))
      .def("get_number_of_cells", &DensityGrid::get_number_of_cells);
}
