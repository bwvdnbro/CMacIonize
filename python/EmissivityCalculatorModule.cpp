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
 * @file EmissivityCalculatorModule.cpp
 *
 * @brief Python module exposure of EmissivityCalculator.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "EmissivityCalculator.hpp"
#include "DensityGrid.hpp"
#include "LineCoolingData.hpp"
#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/module.hpp>
#include <boost/python/numeric.hpp>

/*! @brief Tell numpy to use the non deprecated API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

/**
 * @brief Python constructor for an EmissivityCalculator.
 *
 * @param abundances Abundances of all elements other than hydrogen.
 * @return boost::shared_ptr to the EmissivityCalculator.
 */
static boost::shared_ptr< EmissivityCalculator >
initEmissivityCalculator(Abundances &abundances) {
  return boost::shared_ptr< EmissivityCalculator >(
      new EmissivityCalculator(abundances));
}

/**
 * @brief Python wrapper around EmissivityCalculator::calculate_emissivities().
 *
 * @param calculator EmissivityCalculator on which to act (acts as self).
 * @param grid DensityGrid for which the emissivities are calculated.
 * @return Python dict containing the emissivity values as numpy.ndarrays.
 */
static boost::python::dict get_emissivities(EmissivityCalculator &calculator,
                                            DensityGrid &grid) {
  npy_intp size = grid.get_number_of_cells();
  PyObject *narr = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  boost::python::dict result;
  for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
    result[EmissivityValues::get_name(line)] = arr.copy();
  }

  std::vector< EmissivityValues > emissivities =
      calculator.get_emissivities(grid);

  for (size_t i = 0; i < emissivities.size(); ++i) {
    for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
      result[EmissivityValues::get_name(line)][i] =
          emissivities[i].get_emissivity(line);
    }
  }

  return result;
}

/**
 * @brief Compute the emissivities for the given DensityGrid and store them in
 * the grid itself.
 *
 * @param calculator EmissivityCalculator on which to act (acts as self).
 * @param grid DensityGrid on which to act.
 */
static void compute_emissivities(EmissivityCalculator &calculator,
                                 DensityGrid &grid) {
  calculator.calculate_emissivities(grid);
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
  } else if (name == "OI_6300") {
    return "no idea";
  } else if (name == "Temperature") {
    return "K";
  } else {
    cmac_error("Unknown variable: %s!", name.c_str());
    return "";
  }
}

/**
 * @brief Get the EmissionLine corresponding to the given std::string.
 *
 * @param name Name of an EmissionLine.
 * @return Corresponding EmissionLine.
 */
static EmissionLine get_line(std::string name) {
  if (name == "OI_6300") {
    return EMISSIONLINE_OI_6300;
  } else {
    cmac_error("Unknown emission line name: %s!", name.c_str());
    return NUMBER_OF_EMISSIONLINES;
  }
}

/**
 * @brief Make an emission map for the given emission line in the given
 * direction.
 *
 * @param calculator EmissivityCalculator on which to act (acts as self).
 * @param grid DensityGrid on which to act.
 * @param direction Coordinate direction along which the projection is made.
 * @param line EmissionLine to make a map of.
 * @param shape Size of the resulting map.
 * @return Python dict containing a numpy.ndarray with the requested map, and a
 * string representation of the units in which the variables are expressed.
 */
static boost::python::dict make_emission_map(EmissivityCalculator &calculator,
                                             DensityGrid &grid, char direction,
                                             std::string line,
                                             boost::python::tuple shape) {

  // make sure the grid has emissivity values
  calculator.calculate_emissivities(grid);

  npy_intp size[2] = {boost::python::extract< unsigned int >(shape[0]),
                      boost::python::extract< unsigned int >(shape[1])};
  PyObject *narr = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  Box<> box = grid.get_box();

  for (int_fast32_t i = 0; i < size[0]; ++i) {
    for (int_fast32_t j = 0; j < size[1]; ++j) {
      CoordinateVector<> start;
      CoordinateVector<> dir;
      if (direction == 'x') {
        dir[0] = 1.;
        start[0] = box.get_anchor().x();
        start[1] =
            box.get_anchor().y() + (i + 0.5) * box.get_sides().y() / size[0];
        start[2] =
            box.get_anchor().z() + (j + 0.5) * box.get_sides().z() / size[1];
      } else if (direction == 'y') {
        dir[1] = 1.;
        start[0] =
            box.get_anchor().x() + (i + 0.5) * box.get_sides().x() / size[0];
        start[1] = box.get_anchor().y();
        start[2] =
            box.get_anchor().z() + (j + 0.5) * box.get_sides().z() / size[1];
      } else if (direction == 'z') {
        dir[2] = 1.;
        start[0] =
            box.get_anchor().x() + (i + 0.5) * box.get_sides().x() / size[0];
        start[1] =
            box.get_anchor().y() + (j + 0.5) * box.get_sides().y() / size[1];
        start[2] = box.get_anchor().z();
      } else {
        cmac_error("Unknown coordinate direction: %c!", direction);
      }
      arr[i][j] = grid.get_total_emission(start, dir, get_line(line));
    }
  }

  boost::python::dict result;
  result["values"] = arr.copy();
  result["units"] = get_variable_unit(line);
  return result;
}

/**
 * @brief Python constructor for Abundances.
 *
 * @param abundances Python dict containing abundances for all elements other
 * than hydrogen.
 * @return boost::shared_ptr to an Abundances instance.
 */
static boost::shared_ptr< Abundances >
initAbundances(boost::python::dict abundances) {

  double AHe = boost::python::extract< double >(abundances["helium"]);
  double AC = boost::python::extract< double >(abundances["carbon"]);
  double AN = boost::python::extract< double >(abundances["nitrogen"]);
  double AO = boost::python::extract< double >(abundances["oxygen"]);
  double ANe = boost::python::extract< double >(abundances["neon"]);
  double AS = boost::python::extract< double >(abundances["sulphur"]);
  return boost::shared_ptr< Abundances >(
      new Abundances(AHe, AC, AN, AO, ANe, AS));
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE(libemissivitycalculator) {
  // we need to tell Boost we mean numpy.ndarray whenever we write
  // boost::python::numeric::array
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  // we have to kindly ask numpy to initialize its array functionality
  import_array();

  // we tell Boost we want to expose our version of
  // LineCoolingData.get_line_strengths()
  boost::python::class_< EmissivityCalculator,
                         boost::shared_ptr< EmissivityCalculator > >(
      "EmissivityCalculator", boost::python::no_init)
      .def("__init__",
           boost::python::make_constructor(&initEmissivityCalculator))
      .def("get_emissivities", &get_emissivities)
      .def("compute_emissivities", &compute_emissivities)
      .def("make_emission_map", &make_emission_map);

  boost::python::class_< Abundances, boost::shared_ptr< Abundances >,
                         boost::noncopyable >("Abundances",
                                              boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&initAbundances));
}
