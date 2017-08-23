#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file gamma_NeII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Griffin, D. C, Mitnik, D. M., Badnell, N. R. 2001, JPhB, 34, 4401
# (http://adsabs.harvard.edu/abs/2001JPhB...34.4401G).
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

## load modules
import numpy as np
# for curve_fit
import scipy.optimize as opt
# for plotting (using a backend that does not require a graphics environment)
import matplotlib
matplotlib.use("Agg")
import pylab as pl
# for the fitting curve
from fitting_curve import fitting_curve, print_fit_variables, \
                          initialize_data_values, append_data_values, \
                          print_data_values, get_code, jacobian_fitting_curve, \
                          round_parameters

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Griffin, Mitnik & Badnell (2001), table 4
  T = np.array([1000., 4000., 10000., 40000., 100000., 400000.])
  data = np.array([0.266, 0.299, 0.314, 0.35, 0.4, 0.473])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  imin = 0
  imax = len(T)
  # fit the curve
  A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[imin:imax],
                      maxfev = 1000000, method = "trf",
                      p0 = (0., 1., 1., 1., 0.01, 0.01, 0.),
                      jac = jacobian_fitting_curve)
  A = round_parameters(*A)
  # compute the xi2 difference between the data values (in the fitting
  # interval) and the curve
  xi2 = sum( (data[imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
  # output some info
  print "Transition: 0 to 1"
  print_fit_variables(*A)
  print "convergence:", xi2
  print "validity: [", T[imin], ",", T[imax-1], "]"
  # write the fitting code for this transition
  code += get_code("NeII", "REMOVE_THIS_BLOCK", *A)
  # add the values to the list strings
  append_data_values(data_values, *A)

  # plot the data and fit for visual comparison
  Trange = np.logspace(3., 5., 100)
  pl.plot(T, data, "k.")
  pl.plot(Trange, fitting_curve(Trange, *A), "r-")
  pl.xlim(0., 1.e5)
  pl.savefig("tmp/NeII_{key}.png".format(key = "G0t1"))
  pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
