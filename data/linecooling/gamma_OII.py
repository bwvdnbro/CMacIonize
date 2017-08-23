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
# @file gamma_OII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Kisielius, R., Storey, P. J., Ferland, G. J. & Keenan, F. P. 2009,
# MNRAS, 397, 903 (http://adsabs.harvard.edu/abs/2009MNRAS.397..903K).
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

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Kisielius et al. (2009), table 3
  T = np.array([100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 3000.,
                5000., 7500., 10000., 15000., 20000., 30000., 50000., 75000.,
                100000.])
  # 4S3/2 to 2D5/2
  data["G0t1"] = np.array([0.796, 0.797, 0.798, 0.801, 0.808, 0.817, 0.823,
                           0.830, 0.832, 0.832, 0.831, 0.833, 0.834, 0.839,
                           0.844, 0.856, 0.881, 0.905, 0.919])
  # 4S3/2 to 2D3/2
  data["G0t2"] = np.array([0.531, 0.533, 0.533, 0.535, 0.540, 0.546, 0.550,
                           0.554, 0.555, 0.554, 0.553, 0.553, 0.554, 0.557,
                           0.561, 0.569, 0.585, 0.601, 0.611])
  # 2D5/2 to 2D3/2
  data["G1t2"] = np.array([1.095, 1.086, 1.078, 1.072, 1.097, 1.151, 1.194,
                           1.239, 1.254, 1.256, 1.241, 1.221, 1.203, 1.183,
                           1.179, 1.193, 1.229, 1.257, 1.270])
  # 2P3/2 to 2P1/2
  data["G3t4"] = np.array([0.273, 0.274, 0.274, 0.274, 0.274, 0.275, 0.275,
                           0.276, 0.276, 0.277, 0.279, 0.282, 0.285, 0.294,
                           0.305, 0.327, 0.361, 0.388, 0.405])
  # 4S3/2 to 2P3/2
  data["G0t3"] = np.array([0.244, 0.245, 0.245, 0.245, 0.245, 0.246, 0.246,
                           0.247, 0.247, 0.249, 0.251, 0.253, 0.256, 0.260,
                           0.265, 0.274, 0.290, 0.304, 0.312])
  # 4S3/2 to 2P1/2
  data["G0t4"] = np.array([0.126, 0.126, 0.126, 0.126, 0.127, 0.127, 0.127,
                           0.127, 0.128, 0.128, 0.129, 0.131, 0.132, 0.134,
                           0.136, 0.141, 0.149, 0.155, 0.159])
  # 2D5/2 to 2P3/2
  data["G1t3"] = np.array([0.791, 0.793, 0.793, 0.794, 0.796, 0.797, 0.799,
                           0.801, 0.804, 0.809, 0.820, 0.834, 0.851, 0.891,
                           0.930, 0.997, 1.084, 1.144, 1.178])
  # 2D5/2 to 2P1/2
  data["G1t4"] = np.array([0.315, 0.316, 0.316, 0.316, 0.317, 0.318, 0.318,
                           0.319, 0.320, 0.322, 0.326, 0.332, 0.339, 0.356,
                           0.371, 0.396, 0.427, 0.447, 0.458])
  # 2D3/2 to 2P3/2
  data["G2t3"] = np.array([0.439, 0.440, 0.440, 0.440, 0.441, 0.442, 0.443,
                           0.444, 0.445, 0.448, 0.454, 0.462, 0.472, 0.494,
                           0.516, 0.551, 0.595, 0.624, 0.639])
  # 2D3/2 to 2P1/2
  data["G2t4"] = np.array([0.308, 0.308, 0.309, 0.309, 0.310, 0.310, 0.311,
                           0.312, 0.313, 0.315, 0.319, 0.324, 0.331, 0.345,
                           0.360, 0.386, 0.421, 0.445, 0.459])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    imin = 9
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 1000000,
                        p0 = (0., 1., 100., 1., 0.01, 0.01, 0.),
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
    code += get_code("OII", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.xlim(0., 1.e5)
    pl.savefig("tmp/OII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
