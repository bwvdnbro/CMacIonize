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
#include "DensityGrid.hpp"
#include "EmissivityCalculator.hpp"
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
  for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
    EmissionLine line = static_cast< EmissionLine >(i);
    result[EmissivityValues::get_name(line)] = arr.copy();
  }

  std::vector< EmissivityValues > emissivities =
      calculator.get_emissivities(grid);

  for (unsigned int i = 0; i < emissivities.size(); ++i) {
    for (int j = 0; j < NUMBER_OF_EMISSIONLINES; ++j) {
      EmissionLine line = static_cast< EmissionLine >(j);
      result[EmissivityValues::get_name(line)][i] =
          emissivities[i].get_emissivity(line);
    }
  }

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
  double AS = boost::python::extract< double >(abundances["sulfur"]);
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

  // we tell Boost we want to expose our version of LineCoolingData.linestr()
  boost::python::class_< EmissivityCalculator,
                         boost::shared_ptr< EmissivityCalculator > >(
      "EmissivityCalculator", boost::python::no_init)
      .def("__init__",
           boost::python::make_constructor(&initEmissivityCalculator))
      .def("get_emissivities", &get_emissivities);

  boost::python::class_< Abundances, boost::shared_ptr< Abundances >,
                         boost::noncopyable >("Abundances",
                                              boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&initAbundances));
}
