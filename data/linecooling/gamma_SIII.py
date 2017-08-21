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
# @file gamma_SIII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Hudson, C. E., Ramsbottom, C. A. & Scott, M. P. 2012, ApJ, 750, 65
# (http://adsabs.harvard.edu/abs/2012ApJ...750...65H).
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
                          print_data_values, get_code

# dictionary that links abbreviated transition names to the full names used in
# LineCoolingData
transitions = {
  "G0t1": "TRANSITION_0_to_1",
  "G0t2": "TRANSITION_0_to_2",
  "G0t3": "TRANSITION_0_to_3",
  "G0t4": "TRANSITION_0_to_4",
  "G1t2": "TRANSITION_1_to_2",
  "G1t3": "TRANSITION_1_to_3",
  "G1t4": "TRANSITION_1_to_4",
  "G2t3": "TRANSITION_2_to_3",
  "G2t4": "TRANSITION_2_to_4",
  "G3t4": "TRANSITION_3_to_4"
}

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Hudson, Ramsbottom & Scott (2012), table 5
  logT = np.array([3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2,
                   5.4, 5.6, 5.8, 6.])
  T = 10.**logT
  # 3P0 to 3P1
  data["G0t1"] = np.array([2.08, 2.10, 2.13, 2.20, 2.27, 2.26, 2.17, 2.07, 1.97,
                           1.84, 1.63, 1.36, 1.07, 0.81, 0.59, 0.41])
  # 3P0 to 3P2
  data["G0t2"] = np.array([9.74e-1, 9.49e-1, 9.36e-1, 9.50e-1, 9.73e-1, 1.02,
                           1.11, 1.23, 1.31, 1.31, 1.22, 1.06, 8.74e-1, 6.96e-1,
                           5.45e-1, 4.23e-1])
  # 3P1 to 3P2
  data["G1t2"] = np.array([4.82, 4.80, 4.79, 4.90, 5.03, 5.10, 5.19, 5.31, 5.36,
                           5.19, 4.73, 4.03, 3.25, 2.52, 1.92, 1.44])
  # 1D2 to 1S0
  data["G3t4"] = np.array([8.63e-1, 8.66e-1, 8.72e-1, 9.50e-1, 1.14, 1.38, 1.60,
                           1.79, 1.94, 1.98, 1.90, 1.73, 1.54, 1.35, 1.18,
                           1.02])
  # 3P0 to 1D2
  data["G0t3"] = np.array([6.98e-1, 7.33e-1, 7.38e-1, 7.20e-1, 7.10e-1, 7.29e-1,
                           7.65e-1, 7.91e-1, 7.86e-1, 7.44e-1, 6.60e-1, 5.48e-1,
                           4.30e-1, 3.23e-1, 2.34e-1, 1.64e-1])
  # 3P0 to 1S0
  data["G0t4"] = np.array([8.27e-2, 8.89e-2, 9.58e-2, 1.04e-1, 1.14e-1, 1.25e-1,
                           1.38e-1, 1.54e-1, 1.71e-1, 1.78e-1, 1.72e-1, 1.55e-1,
                           1.32e-1, 1.09e-1, 8.72e-2, 6.84e-2])
  # 3P1 to 1D2
  data["G1t3"] = np.array([2.09, 2.19, 2.19, 2.14, 2.11, 2.17, 2.28, 2.36, 2.35,
                           2.23, 1.98, 1.64, 1.29, 9.65e-1, 6.98e-1, 4.90e-1])
  # 3P1 to 1S0
  data["G1t4"] = np.array([2.17e-1, 2.36e-1, 2.57e-1, 2.79e-1, 3.05e-1, 3.30e-1,
                           3.54e-1, 3.84e-1, 4.08e-1, 4.05e-1, 3.65e-1, 3.01e-1,
                           2.31e-1, 1.69e-1, 1.18e-1, 8.04e-2])
  # 3P2 to 1D2
  data["G2t3"] = np.array([3.90, 4.07, 4.08, 3.99, 3.94, 4.02, 4.19, 4.30, 4.26,
                           4.01, 3.55, 2.95, 2.32, 1.75, 1.27, 8.99e-1])
  # 3P2 to 1S0
  data["G2t4"] = np.array([3.54e-1, 3.87e-1, 4.21e-1, 4.60e-1, 5.04e-1, 5.45e-1,
                           5.85e-1, 6.34e-1, 6.73e-1, 6.66e-1, 5.99e-1, 4.94e-1,
                           3.80e-1, 2.77e-1, 1.95e-1, 1.33e-1])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    # we force the curve to go through the value at 10,000 K using a global
    # variable
    norm = data[key][5]
    # we start by fitting to the full data set
    imin = 0
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 1000000,
                        p0 = (1., 1., 100., 1., 1., 1., 1., 1.))
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum( (data[key][imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
    # output some info
    print "Transition:", key
    print_fit_variables(*A)
    print "convergence:", xi2
    print "validity: [", T[imin], ",", T[imax-1], "]"
    # write the fitting code for this transition
    code += get_code("CIII", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.xlim(0., 1.e5)
    pl.savefig("tmp/SIII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
