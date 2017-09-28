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
# @file gamma_OIII.py
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
                          print_data_values, get_code, jacobian_fitting_curve, \
                          round_parameters

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
  data["G0t1"] = np.array([0.4975, 0.5066, 0.5115, 0.518, 0.5296, 0.5454,
                           0.559, 0.5678, 0.5788, 0.5918, 0.5938])
  # 3P0 to 3P2
  data["G0t2"] = np.array([0.2455, 0.2493, 0.2509, 0.2541, 0.2609, 0.2713,
                           0.2832, 0.2955, 0.3101, 0.3254, 0.3314])
  # 3P1 to 3P2
  data["G1t2"] = np.array([1.173, 1.193, 1.203, 1.218, 1.248, 1.291,
                           1.335, 1.373, 1.419, 1.468, 1.482])
  # 3P? to 1D2
  GPtD = np.array([2.2233, 2.1888, 2.1416, 2.1117, 2.1578, 2.2892,
                   2.4497, 2.5851, 2.673, 2.7019, 2.6594])
  # 3P? to 1S0
  GPtS = np.array([0.2754, 0.2738, 0.2713, 0.2693, 0.2747, 0.2925,
                   0.3174, 0.3405, 0.3563, 0.3621, 0.3571])
  # 1D2 to 1S0
  data["G3t4"] = np.array([0.4241, 0.4268, 0.4357, 0.4652, 0.5232, 0.5815,
                           0.61, 0.609, 0.5971, 0.5865, 0.5725])

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
    imin = 2
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 1000000,
                        p0 = (0., 1., 1., 1., 0.1, 0.1, 0.),
                        jac = jacobian_fitting_curve)
    A = round_parameters(*A)
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum( (data[key][imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
    # output some info
    print "Transition:", key
    print_fit_variables(*A)
    print "convergence:", xi2
    print "validity: [", T[imin], ",", T[imax-1], "]"
    # write the fitting code for this transition
    code += get_code("OIII", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.savefig("tmp/OIII_{key}.png".format(key = key))
    pl.close()

    # save the plot values in separate files
    dfile = open("tmp/OIII_{key}_data.txt".format(key = key), "w")
    for i in range(len(T)):
      dfile.write("{T}\t{data}\n".format(T = T[i], data = data[key][i]))
    dfile.close()
    ffile = open("tmp/OIII_{key}_fit.txt".format(key = key), "w")
    for i in range(len(Trange)):
      ffile.write("{T}\t{fit}\n".format(T = Trange[i],
                                        fit = fitting_curve(Trange[i], *A)))
    ffile.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
