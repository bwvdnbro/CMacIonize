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
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/module.hpp>
#include <cmath>

boost::python::list python_linestr(LineCoolingData &lines,
                                   boost::python::list abundances) {
  double abund[12];
  for (unsigned int i = 0; i < 12; ++i) {
    abund[i] = boost::python::extract< double >(abundances[i]);
  }

  double c6300 = 0., c9405 = 0., c6312 = 0., c33mu = 0., c19mu = 0., c3729 = 0.,
         c3727 = 0., c7330 = 0., c4363 = 0., c5007 = 0., c52mu = 0., c88mu = 0.,
         c5755 = 0., c6584 = 0., c4072 = 0., c6717 = 0., c6725 = 0., c3869 = 0.,
         cniii57 = 0., cneii12 = 0., cneiii15 = 0., cnii122 = 0., cii2325 = 0.,
         ciii1908 = 0., coii7325 = 0., csiv10 = 0.;
  lines.linestr(8000., 100., abund, c6300, c9405, c6312, c33mu, c19mu, c3729,
                c3727, c7330, c4363, c5007, c52mu, c88mu, c5755, c6584, c4072,
                c6717, c6725, c3869, cniii57, cneii12, cneiii15, cnii122,
                cii2325, ciii1908, coii7325, csiv10);

  boost::python::list linestrengths;
  linestrengths.append(c6300);
  linestrengths.append(c9405);
  linestrengths.append(c6312);
  linestrengths.append(c33mu);
  linestrengths.append(c19mu);
  linestrengths.append(c3729);
  linestrengths.append(c3727);
  linestrengths.append(c7330);
  linestrengths.append(c4363);
  linestrengths.append(c5007);
  linestrengths.append(c52mu);
  linestrengths.append(c88mu);
  linestrengths.append(c5755);
  linestrengths.append(c6584);
  linestrengths.append(c4072);
  linestrengths.append(c6717);
  linestrengths.append(c6725);
  linestrengths.append(c3869);
  linestrengths.append(cniii57);
  linestrengths.append(cneii12);
  linestrengths.append(cneiii15);
  linestrengths.append(cnii122);
  linestrengths.append(cii2325);
  linestrengths.append(ciii1908);
  linestrengths.append(coii7325);
  linestrengths.append(csiv10);

  return linestrengths;
}

BOOST_PYTHON_MODULE(liblinecoolingdata) {
  boost::python::class_< LineCoolingData >("LineCoolingData")
      .def("linestr", &python_linestr);
}
