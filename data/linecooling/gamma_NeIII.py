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
# @file gamma_NeIII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Butler, K. & Zeippen, C. J. 1994, A&AS, 108, 1
# (http://adsabs.harvard.edu/abs/1994A%26AS..108....1B).
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
initial_guess = (1., 1., 100., 1., 1., 0.1, 1.)
transitions = {
  "G0t1": ["TRANSITION_0_to_1",
           (0.348500741382, 0.195852073826, -36.4501581617, 0.,
            -0.0181252763926, 0., 1.00000001306)],
  "G0t2": ["TRANSITION_0_to_2",
           (0.242259790973, 0.118399224056, -28.6257166059, 0.,
            -0.010466296122, 0., 1.00000000994)],
  "G0t3": ["TRANSITION_0_to_3",
           (0.994546627597, -9.74092739212e-06, 0.813424494789, 0.,
            7.92513938776e-07, 0., 1.00000000001)],
  "G0t4": ["TRANSITION_0_to_4",
           (1.10193063883, -8.06658428298e-06, 0.0418955574988, 0.,
            8.07224313295e-07, 0., 0.999999999982)],
  "G1t2": ["TRANSITION_1_to_2",
           (0.382673053848, 0.0454578976384, -8.08652985413, 0.,
            -0.00421441235504, 0., 0.99999999998)],
  "G1t3": ["TRANSITION_1_to_3",
           (0.998072239297, -8.15868142461e-06, 0.475815954012, 0.,
            7.01067249995e-07, 0., 0.999999999989)],
  "G1t4": ["TRANSITION_1_to_4",
           (0.192021167944, 0.0357591526476, -1.80494222891, 0.,
            -0.00313573333494, 0., 0.999999998745)],
  "G2t3": ["TRANSITION_2_to_3",
           (0.9788720275, 7.00840864284e-06, 0.178725598354, 0.,
            -7.70012994119e-07, 0., 1.0)],
  "G2t4": ["TRANSITION_2_to_4",
           (0.370195520847, 0.00309656666641, 0.209609657805, 0.,
            -0.000292624490616, 0., 1.00000000001)],
  "G3t4": ["TRANSITION_3_to_4",
           (1.05817326234, -4.05953358177e-05, 0.187708099593, 0.,
            4.31174660385e-06, 0., 1.0)]
}

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Butler & Zeippen (1994), table 2
  logT = np.array([3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.])
  T = 10.**logT
  # 3P2 to 3P1
  data["G0t1"] = np.array([0.481, 0.545, 0.634, 0.708, 0.752, 0.774, 0.778,
                           0.771, 0.764, 0.773, 0.794])
  # 3P2 to 3P0
  data["G0t2"] = np.array([0.128, 0.149, 0.174, 0.194, 0.204, 0.208, 0.208,
                           0.205, 0.203, 0.205, 0.211])
  # 3P1 to 3P0
  data["G1t2"] = np.array([0.154, 0.168, 0.194, 0.218, 0.235, 0.244, 0.247,
                           0.247, 0.246, 0.249, 0.256])
  # 1D2 to 1S0
  data["G3t4"] = np.array([0.266, 0.266, 0.266, 0.266, 0.267, 0.269, 0.277,
                           0.292, 0.31, 0.325, 0.333])
  # 3P2 to 1D2
  data["G0t3"] = np.array([0.749, 0.765, 0.771, 0.767, 0.76, 0.754, 0.749,
                           0.746, 0.749, 0.763, 0.782])
  # 3P2 to 1S0
  data["G0t4"] = np.array([0.083, 0.083, 0.083, 0.083, 0.084, 0.084, 0.084,
                           0.086, 0.09, 0.095, 0.1])
  # 3P1 to 1D2
  data["G1t3"] = np.array([0.45, 0.459, 0.462, 0.46, 0.456, 0.452, 0.449,
                           0.448, 0.449, 0.458, 0.469])
  # 3P1 to 1S0
  data["G1t4"] = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.051,
                           0.054, 0.057, 0.06])
  # 3P0 to 1D2
  data["G2t3"] = np.array([0.15, 0.153, 0.154, 0.153, 0.152, 0.151, 0.15,
                           0.149, 0.15, 0.153, 0.156])
  # 3P0 to 1S0
  data["G2t4"] = np.array([0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017,
                           0.017, 0.018, 0.019, 0.02])

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
                        p0 = transitions[key][1],
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
    code += get_code("NeIII", transitions[key][0], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.savefig("tmp/NeIII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
