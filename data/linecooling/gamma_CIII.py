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
# @file gamma_CIII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Berrington, K. A., Burke, P. G., Dufton, P. L. & Kingston, A. E.
# 1985, ADNDT, 33, 195 (http://adsabs.harvard.edu/abs/1985ADNDT..33..195B).
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

  # data from Berrington et al. (1985), table X
  logT = np.array([4.1, 4.3, 4.5, 4.7, 4.9, 5.1, 5.3, 5.5, 5.7])
  T = 10.**logT

  # 1S0 to 3P?
  GStP = np.array([1.04, 1.03, 1.03, 9.96e-1, 9.16e-1, 8.08e-1, 6.92e-1,
                   5.79e-1, 4.73e-1])
  # 3P? to 1P1
  GPtP = np.array([3.7, 3.54, 3.29, 2.93, 2.5, 2.07, 1.68, 1.38, 1.2])
  # 1S0 to 3P0
  data["G0t1"] = GStP / 9.
  # 1S0 to 3P1
  data["G0t2"] = GStP / 3.
  # 3P0 to 3P1
  data["G1t2"] = np.array([9.62e-1, 1.03, 1.09, 1.11, 1.06, 9.53e-1, 8.15e-1,
                           6.73e-1, 5.41e-1])
  # 3P2 to 1P1
  data["G3t4"] = 5.* GPtP / 9.
  # 1S0 to 3P2
  data["G0t3"] = 5. * GStP / 9.
  # 1S0 to 1P1
  data["G0t4"] = np.array([4., 4.11, 4.24, 4.42, 4.65, 5., 5.5, 6.17, 7.01])
  # 3P0 to 3P2
  data["G1t3"] = np.array([7.18e-1, 8.57e-1, 1.01, 1.11, 1.09, 9.88e-1, 8.43e-1,
                           7.03e-1, 5.88e-1])
  # 3P0 to 1P1
  data["G1t4"] = GPtP / 9.
  # 3P1 to 3P2
  data["G2t3"] = np.array([2.78, 3.17, 3.6, 3.83, 3.74, 3.37, 2.88, 2.39, 1.97])
  # 3P1 to 1P1
  data["G2t4"] = GPtP / 3.

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    # not exactly right, but good enough
    norm = data[key][0]
    # we start by fitting to the full data set
    imin = 0
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax])
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum( (data[key][imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
    # if xi2 is too large: shrink the fitting interval and try again
    while xi2 > 1.e100:
      imin += 1
      imax -= 1
      # if the interval becomes too small, we bail out
      if imax - imin < 2:
        break
      A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax])
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
    pl.savefig("tmp/CIII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
