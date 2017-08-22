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
# @file gamma_NII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Lennon, D. J. & Burke, V. M. 1994, A&AS, 103, 273
# (http://adsabs.harvard.edu/abs/1994A%26AS..103..273L).
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
                          print_data_values, get_code, jacobian_fitting_curve

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

# main function: computes fits to the data of Lennon & Burke (1994) and plots
# the data and fits for visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Lennon & Burke (1994), table 2
  logT = np.array([3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.])
  T = 10.**logT
  # 3P0 to 3P1
  data["G0t1"] = np.array([0.3336, 0.3333, 0.3394, 0.3581, 0.3837, 0.4076,
                           0.4293, 0.4566, 0.491, 0.5238, 0.5499])
  # 3P0 to 3P2
  data["G0t2"] = np.array([0.2142, 0.2194, 0.2256, 0.2351, 0.2499, 0.272,
                           0.3012, 0.3316, 0.3542, 0.3572, 0.3494])
  # 3P1 to 3P2
  data["G1t2"] = np.array([0.8978, 0.9092, 0.9307, 0.9753, 1.041, 1.12, 1.213,
                           1.316,  1.405, 1.458, 1.473])
  # 3P? to 1D2
  GPtD = np.array([2.4841, 2.4963, 2.5163, 2.5469, 2.5884, 2.6409, 2.7,
                   2.7663, 2.8444, 2.9235, 2.9668])
  # 3P? to 1S0
  GPtS = np.array([0.2813, 0.2821, 0.2833, 0.2853, 0.2884, 0.2931,
                   0.2999, 0.3099, 0.324, 0.3396, 0.3524])
  # 1D2 to 1S0
  data["G3t4"] = np.array([1.164, 1.126, 1.071, 0.9998, 0.9176, 0.8338, 0.7605,
                           0.7097, 0.6825, 0.6686, 0.6571])

  # formula (1) in Lennon & Burke (1994)
  data["G0t3"] = GPtD/9.
  data["G0t4"] = GPtS/9.
  data["G1t3"] = GPtD/3.
  data["G1t4"] = GPtS/3.
  data["G2t3"] = 5.*GPtD/9.
  data["G2t4"] = 5.*GPtS/9.

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    imin = 0
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 100000,
                        p0 = (1., 1., 10., 0., 1., 0., 1.),
                        jac = jacobian_fitting_curve)
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum( (data[key][imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
    # output some info
    print "Transition:", key
    print_fit_variables(*A)
    print "convergence:", xi2
    print "validity: [", T[imin], ",", T[imax-1], "]"
    # write the fitting code for this transition
    code += get_code("NII", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.savefig("tmp/NII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
