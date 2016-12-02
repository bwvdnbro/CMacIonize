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
 * @file LineCoolingDataModule.cpp
 *
 * @brief Python module exposure of LineCoolingData.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Error.hpp"
#include "LineCoolingData.hpp"
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/lvalue_from_pytype.hpp>
#include <boost/python/module.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/object.hpp>
#include <cmath>

/*! @brief Tell numpy to use the non deprecated API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

/**
 * @brief Python version of LineCoolingData::linestr().
 *
 * @param lines LineCoolingData object that is wrapped by this function.
 * @param T Temperatures (in K).
 * @param ne Electron densities (in m^-3).
 * @param abundances boost::python::list of abundances.
 * @return boost::python::list containing the output of
 * LineCoolingData::linestr().
 */
static boost::python::object python_linestr(LineCoolingData &lines,
                                            boost::python::numeric::array &T,
                                            boost::python::numeric::array &ne,
                                            boost::python::list abundances) {
  double abund[12];
  for (unsigned int i = 0; i < 12; ++i) {
    abund[i] = boost::python::extract< double >(abundances[i]);
  }

  boost::python::tuple Tshape =
      boost::python::extract< boost::python::tuple >(T.attr("shape"));
  const unsigned int numT = boost::python::len(Tshape);

  npy_intp size[2] = {numT, 26};
  PyObject *narr = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);

  for (unsigned int iT = 0; iT < numT; ++iT) {
    double Ti = boost::python::extract< double >(T[iT]);
    double nei = boost::python::extract< double >(ne[iT]);

    double c6300 = 0., c9405 = 0., c6312 = 0., c33mu = 0., c19mu = 0.,
           c3729 = 0., c3727 = 0., c7330 = 0., c4363 = 0., c5007 = 0.,
           c52mu = 0., c88mu = 0., c5755 = 0., c6584 = 0., c4072 = 0.,
           c6717 = 0., c6725 = 0., c3869 = 0., cniii57 = 0., cneii12 = 0.,
           cneiii15 = 0., cnii122 = 0., cii2325 = 0., ciii1908 = 0.,
           coii7325 = 0., csiv10 = 0.;
    lines.linestr(Ti, nei, abund, c6300, c9405, c6312, c33mu, c19mu, c3729,
                  c3727, c7330, c4363, c5007, c52mu, c88mu, c5755, c6584, c4072,
                  c6717, c6725, c3869, cniii57, cneii12, cneiii15, cnii122,
                  cii2325, ciii1908, coii7325, csiv10);

    arr[boost::python::make_tuple(iT, 0)] = c6300;
    arr[boost::python::make_tuple(iT, 1)] = c9405;
    arr[boost::python::make_tuple(iT, 2)] = c6312;
    arr[boost::python::make_tuple(iT, 3)] = c33mu;
    arr[boost::python::make_tuple(iT, 4)] = c19mu;
    arr[boost::python::make_tuple(iT, 5)] = c3729;
    arr[boost::python::make_tuple(iT, 6)] = c3727;
    arr[boost::python::make_tuple(iT, 7)] = c7330;
    arr[boost::python::make_tuple(iT, 8)] = c4363;
    arr[boost::python::make_tuple(iT, 9)] = c5007;
    arr[boost::python::make_tuple(iT, 10)] = c52mu;
    arr[boost::python::make_tuple(iT, 11)] = c88mu;
    arr[boost::python::make_tuple(iT, 12)] = c5755;
    arr[boost::python::make_tuple(iT, 13)] = c6584;
    arr[boost::python::make_tuple(iT, 14)] = c4072;
    arr[boost::python::make_tuple(iT, 15)] = c6717;
    arr[boost::python::make_tuple(iT, 16)] = c6725;
    arr[boost::python::make_tuple(iT, 17)] = c3869;
    arr[boost::python::make_tuple(iT, 18)] = cniii57;
    arr[boost::python::make_tuple(iT, 19)] = cneii12;
    arr[boost::python::make_tuple(iT, 20)] = cneiii15;
    arr[boost::python::make_tuple(iT, 21)] = cnii122;
    arr[boost::python::make_tuple(iT, 22)] = cii2325;
    arr[boost::python::make_tuple(iT, 23)] = ciii1908;
    arr[boost::python::make_tuple(iT, 24)] = coii7325;
    arr[boost::python::make_tuple(iT, 25)] = csiv10;
  }

  return arr.copy();
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE(liblinecoolingdata) {
  // we need to tell Boost we mean numpy.ndarray whenever we write
  // boost::python::numeric::array
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  // we have to kindly ask numpy to initialize its array functionality
  import_array();

  // we tell Boost we want to expose our version of LineCoolingData.linestr()
  boost::python::class_< LineCoolingData >("LineCoolingData")
      .def("linestr", &python_linestr);
}
