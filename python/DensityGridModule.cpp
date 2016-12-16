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
#include <boost/python/dict.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/module.hpp>
#include <boost/python/numeric.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

/*! @brief Tell numpy to use the non deprecated API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

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

  // make sure AMR grids are processed correctly
  if (parameters.get_value< std::string >("densitygrid.type") == "AMR") {
    // this overrides whatever value was in that field
    parameters.add_value("densitygrid.amrrefinementscheme.type", "CMacIonize");
  }

  CMacIonizeSnapshotDensityFunction density_function(filename);

  return boost::shared_ptr< DensityGrid >(
      DensityGridFactory::generate(parameters, density_function));
}

/**
 * @brief Get a numpy.ndarray containing the dimensions of the box containing
 * the grid.
 *
 * @param grid DensityGrid on which to act (acts as self).
 * @return Python dict containing two separate arrays for the origin and sides
 * of the box, and a string representation of the units in which the box size is
 * expressed (m).
 */
static boost::python::dict get_box(DensityGrid &grid) {
  Box box = grid.get_box();

  npy_intp size = 3;
  PyObject *narr = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  boost::python::dict result;
  result["origin"] = arr.copy();
  result["sides"] = arr.copy();

  result["origin"][0] = box.get_anchor().x();
  result["origin"][1] = box.get_anchor().y();
  result["origin"][2] = box.get_anchor().z();
  result["sides"][0] = box.get_sides().x();
  result["sides"][1] = box.get_sides().y();
  result["sides"][2] = box.get_sides().z();
  result["units"] = "m";

  return result;
}

/**
 * @brief Get the variable with the given name in the given cell.
 *
 * @param cell DensityValues of a cell.
 * @param name std::string representation of a cell variable name.
 * @return Value of that variable (in SI units).
 */
static double get_single_variable(DensityValues &cell, std::string name) {
  // names are ordered alphabetically
  if (name.find("NeutralFraction") == 0) {
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      if (name == "NeutralFraction" + get_ion_name(i)) {
        IonName ion = static_cast< IonName >(i);
        return cell.get_ionic_fraction(ion);
      }
    }
    cmac_error("Unknown variable: %s!", name.c_str());
    return 0.;
  } else if (name == "NumberDensity") {
    return cell.get_total_density();
  } else if (name == "Temperature") {
    return cell.get_temperature();
  } else {
    cmac_error("Unknown variable: %s!", name.c_str());
    return 0.;
  }
}

/**
 * @brief Get the units in which the given variable name is expressed.
 *
 * @param name std::string representation of a cell variable name.
 * @return std::string representation of the units in which that variable is
 * expressed.
 */
static std::string get_variable_unit(std::string name) {
  // names are ordered alphabetically
  if (name.find("NeutralFraction") == 0) {
    // all neutral fractions are dimensionless
    return "";
  } else if (name == "NumberDensity") {
    return "m^-3";
  } else if (name == "Temperature") {
    return "K";
  } else {
    cmac_error("Unknown variable: %s!", name.c_str());
    return "";
  }
}

/**
 * @brief Get a numpy.ndarray containing the values of the variable with the
 * given name for all cells.
 *
 * @param grid DensityGrid on which to act (acts as self).
 * @param name std::string representation of a cell variable name supported by
 * get_single_variable().
 * @return Python dict containing a numpy.ndarray with the values of the
 * variable for all cells, and a string representation of the units in which the
 * variable is expressed.
 */
static boost::python::dict get_variable(DensityGrid &grid, std::string name) {
  npy_intp size = grid.get_number_of_cells();
  PyObject *narr = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  unsigned int index = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    arr[index] = get_single_variable(it.get_values(), name);
    ++index;
  }

  boost::python::dict result;
  result["values"] = arr.copy();
  result["units"] = get_variable_unit(name);

  return result;
}

/**
 * @brief Get a numpy.ndarray containing the values of the variable with the
 * given name for a planar cut parallel to two of the three coordinate axes.
 *
 * @param grid DensityGrid on which to act (acts as self).
 * @param name std::string representation of a cell variable name supported by
 * get_single_variable().
 * @param coordinate Coordinate axis not parallel to the cut plane (possible
 * values: x, y, or z).
 * @param intercept Value of the intersection point with the axis and the cut
 * plane.
 * @param shape Shape of the 2D result array, with the values corresponding to
 * the lowest parallel axis being in the rows.
 * @return Python dict containing a numpy.ndarray with the requested values, and
 * a string representation of the units in which the variables are expressed.
 */
static boost::python::dict get_variable_cut(DensityGrid &grid, std::string name,
                                            char coordinate, double intercept,
                                            boost::python::tuple shape) {
  if (coordinate != 'x' && coordinate != 'y' && coordinate != 'z') {
    cmac_error("Unknown coordinate: %c!", coordinate);
  }

  npy_intp size[2] = {boost::python::extract< unsigned int >(shape[0]),
                      boost::python::extract< unsigned int >(shape[1])};
  PyObject *narr = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  Box box = grid.get_box();

  double di, dj;
  if (coordinate == 'x') {
    di = box.get_sides().y() / size[0];
    dj = box.get_sides().z() / size[1];
  } else if (coordinate == 'y') {
    di = box.get_sides().x() / size[0];
    dj = box.get_sides().z() / size[1];
  } else {
    di = box.get_sides().x() / size[0];
    dj = box.get_sides().y() / size[1];
  }
  for (unsigned int i = 0; i < size[0]; ++i) {
    for (unsigned int j = 0; j < size[1]; ++j) {
      CoordinateVector<> position;
      if (coordinate == 'x') {
        position[0] = intercept;
        position[1] = box.get_anchor().y() + (i + 0.5) * di;
        position[2] = box.get_anchor().z() + (j + 0.5) * dj;
      } else if (coordinate == 'y') {
        position[0] = box.get_anchor().x() + (i + 0.5) * di;
        position[1] = intercept;
        position[2] = box.get_anchor().z() + (j + 0.5) * dj;
      } else {
        position[0] = box.get_anchor().x() + (i + 0.5) * di;
        position[1] = box.get_anchor().y() + (j + 0.5) * dj;
        position[2] = intercept;
      }
      DensityValues &cell = grid.get_cell_values(position);
      arr[i][j] = get_single_variable(cell, name);
    }
  }

  boost::python::dict result;
  result["values"] = arr.copy();
  result["units"] = get_variable_unit(name);

  return result;
}

/**
 * @brief Get a numpy.ndarray containing the coordinates of all cells in the
 * grid.
 *
 * @param grid DensityGrid on which to act (acts as self).
 * @return Python dict containing a numpy.ndarray with the values of the
 * coordinates for all cells, and a string representation of the units in which
 * the coordinates are expressed (m).
 */
static boost::python::dict get_coordinates(DensityGrid &grid) {
  npy_intp size[2] = {grid.get_number_of_cells(), 3};
  PyObject *narr = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  unsigned int index = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    CoordinateVector<> coords = it.get_cell_midpoint();
    arr[index][0] = coords.x();
    arr[index][1] = coords.y();
    arr[index][2] = coords.z();
    ++index;
  }

  boost::python::dict result;
  result["values"] = arr.copy();
  result["units"] = "m";

  return result;
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE(libdensitygrid) {
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  import_array();

  // we need to use no_init to tell Boost that we provide a custom constructor.
  // we need to use noncopyable to tell Boost that we want to construct an
  // object of an abstract type.
  boost::python::class_< DensityGrid, boost::shared_ptr< DensityGrid >,
                         boost::noncopyable >("DensityGrid",
                                              boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&initDensityGrid))
      .def("get_number_of_cells", &DensityGrid::get_number_of_cells)
      .def("get_box", &get_box)
      .def("get_variable", &get_variable)
      .def("get_variable_cut", &get_variable_cut)
      .def("get_coordinates", &get_coordinates);
}
