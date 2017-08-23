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
#include "LineCoolingData.hpp"
#include "Error.hpp"
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
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
 * @param T numpy.ndarray containing temperatures (in K).
 * @param ne numpy.ndarray containing electron densities (in m^-3).
 * @param abundances boost::python::list of abundances.
 * @return boost::python::dict containing, for every line strength returned by
 * LineCoolingData.linestr(), a numpy.ndarray containing the values for each
 * temperature and electron density.
 */
static boost::python::dict python_get_line_strengths(
    LineCoolingData &lines, boost::python::numeric::array &T,
    boost::python::numeric::array &ne, boost::python::list abundances) {
  double abund[12];
  for (unsigned int i = 0; i < 12; ++i) {
    abund[i] = boost::python::extract< double >(abundances[i]);
  }

  // retrieve numpy.ndarray shape info and check if it is what we expect
  boost::python::tuple Tshape =
      boost::python::extract< boost::python::tuple >(T.attr("shape"));
  const unsigned int Tdim = boost::python::len(Tshape);
  if (Tdim != 1) {
    cmac_error("Expected a 1D array, but got a %uD array!", Tdim);
  }
  boost::python::tuple neshape =
      boost::python::extract< boost::python::tuple >(ne.attr("shape"));
  const unsigned int nedim = boost::python::len(neshape);
  if (nedim != 1) {
    cmac_error("Expected a 1D-array, but got a %uD array!", nedim);
  }
  const unsigned int numT = boost::python::extract< unsigned int >(Tshape[0]);
  const unsigned int numne = boost::python::extract< unsigned int >(neshape[0]);
  if (numT != numne) {
    cmac_error("Temperature and electron density arrays have different sizes "
               "(len(T) = %u, len(ne) = %u)!",
               numT, numne);
  }

  boost::python::dict result;

  // we create a single ndarray with the correct size
  npy_intp size = numT;
  PyObject *narr = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
  boost::python::handle<> handle(narr);
  boost::python::numeric::array arr(handle);
  // we now simply copy that array into the different dictionary elements
  result["c6300"] = arr.copy();
  result["c9405"] = arr.copy();
  result["c6312"] = arr.copy();
  result["c33mu"] = arr.copy();
  result["c19mu"] = arr.copy();
  result["c3729"] = arr.copy();
  result["c3727"] = arr.copy();
  result["c7330"] = arr.copy();
  result["c4363"] = arr.copy();
  result["c5007"] = arr.copy();
  result["c52mu"] = arr.copy();
  result["c88mu"] = arr.copy();
  result["c5755"] = arr.copy();
  result["c6584"] = arr.copy();
  result["c4072"] = arr.copy();
  result["c6717"] = arr.copy();
  result["c6725"] = arr.copy();
  result["c3869"] = arr.copy();
  result["cniii57"] = arr.copy();
  result["cneii12"] = arr.copy();
  result["cneiii15"] = arr.copy();
  result["cnii122"] = arr.copy();
  result["cii2325"] = arr.copy();
  result["ciii1908"] = arr.copy();
  result["coii7325"] = arr.copy();
  result["csiv10"] = arr.copy();

  for (unsigned int iT = 0; iT < numT; ++iT) {
    double Ti = boost::python::extract< double >(T[iT]);
    double nei = boost::python::extract< double >(ne[iT]);

    double c6300 = 0., c9405 = 0., c6312 = 0., c33mu = 0., c19mu = 0.,
           c3729 = 0., c3727 = 0., c7330 = 0., c4363 = 0., c5007 = 0.,
           c52mu = 0., c88mu = 0., c5755 = 0., c6584 = 0., c4072 = 0.,
           c6717 = 0., c6725 = 0., c3869 = 0., cniii57 = 0., cneii12 = 0.,
           cneiii15 = 0., cnii122 = 0., cii2325 = 0., ciii1908 = 0.,
           coii7325 = 0., csiv10 = 0.;
    std::vector< std::vector< double > > line_strengths =
        lines.get_line_strengths(Ti, nei, abund);

    // NII
    c5755 = line_strengths[NII][TRANSITION_3_to_4];
    c6584 = line_strengths[NII][TRANSITION_2_to_3];
    cnii122 = line_strengths[NII][TRANSITION_1_to_2];

    // OI
    c6300 = line_strengths[OI][TRANSITION_0_to_3] +
            line_strengths[OI][TRANSITION_1_to_3];

    // OII
    c3729 = line_strengths[OII][TRANSITION_0_to_1];
    c3727 = line_strengths[OII][TRANSITION_0_to_1] +
            line_strengths[OII][TRANSITION_0_to_2];
    coii7325 = line_strengths[OII][TRANSITION_1_to_4] +
               line_strengths[OII][TRANSITION_2_to_4] +
               line_strengths[OII][TRANSITION_1_to_3] +
               line_strengths[OII][TRANSITION_2_to_3];

    // OIII
    c4363 = line_strengths[OIII][TRANSITION_3_to_4];
    c5007 = line_strengths[OIII][TRANSITION_2_to_3];
    c52mu = line_strengths[OIII][TRANSITION_1_to_2];
    c88mu = line_strengths[OIII][TRANSITION_0_to_1];

    // NeIII
    c3869 = line_strengths[NeIII][TRANSITION_0_to_3];
    cneiii15 = line_strengths[NeIII][TRANSITION_0_to_1];

    // SII
    c4072 = line_strengths[SII][TRANSITION_0_to_3] +
            line_strengths[SII][TRANSITION_0_to_4];
    c6717 = line_strengths[SII][TRANSITION_0_to_2];
    c6725 = line_strengths[SII][TRANSITION_0_to_1] +
            line_strengths[SII][TRANSITION_0_to_2];

    // SIII
    c9405 = line_strengths[SIII][TRANSITION_1_to_3] +
            line_strengths[SIII][TRANSITION_2_to_3];
    c6312 = line_strengths[SIII][TRANSITION_3_to_4];
    c33mu = line_strengths[SIII][TRANSITION_0_to_1];
    c19mu = line_strengths[SIII][TRANSITION_1_to_2];

    // CII
    cii2325 = line_strengths[CII][TRANSITION_0_to_2] +
              line_strengths[CII][TRANSITION_1_to_2] +
              line_strengths[CII][TRANSITION_0_to_3] +
              line_strengths[CII][TRANSITION_1_to_3] +
              line_strengths[CII][TRANSITION_0_to_4] +
              line_strengths[CII][TRANSITION_1_to_4];

    // CIII
    ciii1908 = line_strengths[CIII][TRANSITION_0_to_1] +
               line_strengths[CIII][TRANSITION_0_to_2] +
               line_strengths[CIII][TRANSITION_0_to_3];

    // NIII
    cniii57 = line_strengths[NIII][0];

    // NeII
    cneii12 = line_strengths[NeII][0];

    // not set!!
    c7330 = 0.;
    csiv10 = 0.;

    result["c6300"][iT] = c6300;
    result["c9405"][iT] = c9405;
    result["c6312"][iT] = c6312;
    result["c33mu"][iT] = c33mu;
    result["c19mu"][iT] = c19mu;
    result["c3729"][iT] = c3729;
    result["c3727"][iT] = c3727;
    result["c7330"][iT] = c7330;
    result["c4363"][iT] = c4363;
    result["c5007"][iT] = c5007;
    result["c52mu"][iT] = c52mu;
    result["c88mu"][iT] = c88mu;
    result["c5755"][iT] = c5755;
    result["c6584"][iT] = c6584;
    result["c4072"][iT] = c4072;
    result["c6717"][iT] = c6717;
    result["c6725"][iT] = c6725;
    result["c3869"][iT] = c3869;
    result["cniii57"][iT] = cniii57;
    result["cneii12"][iT] = cneii12;
    result["cneiii15"][iT] = cneiii15;
    result["cnii122"][iT] = cnii122;
    result["cii2325"][iT] = cii2325;
    result["ciii1908"][iT] = ciii1908;
    result["coii7325"][iT] = coii7325;
    result["csiv10"][iT] = csiv10;
  }

  return result;
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

  // we tell Boost we want to expose our version of
  // LineCoolingData.get_line_strengths()
  boost::python::class_< LineCoolingData >("LineCoolingData")
      .def("get_line_strengths", &python_get_line_strengths);
}
