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
# @file gamma_CII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Tayal, S. S. 2008, A&A, 486, 629
# (http://adsabs.harvard.edu/abs/2008A%26A...486..629T).
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

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  data = {}

  # data from Tayal (2008), table 5
  T = np.array([1000., 2000., 3000., 5000., 10000., 20000., 30000., 50000.,
                100000., 200000., 400000.])
  # 2P1/2 to 2P3/2
  data["G0t1"] = np.array([1.551, 1.579, 1.621, 1.744, 2.028, 2.200, 2.193,
                           2.102, 1.874, 1.586, 1.298])
  # 2P1/2 to 4P1/2
  data["G0t2"] = np.array([0.229, 0.253, 0.259, 0.260, 0.261, 0.262, 0.259,
                           0.246, 0.218, 0.185, 0.150])
  # 2P3/2 to 4P1/2
  data["G1t2"] = np.array([0.126, 0.158, 0.167, 0.173, 0.181, 0.188, 0.188,
                           0.179, 0.158, 0.132, 0.107])
  # 4P3/2 to 4P5/2
  data["G3t4"] = np.array([1.565, 1.575, 1.592, 1.658, 1.926, 2.300, 2.439,
                           2.496, 2.443, 2.328, 2.146])
  # 2P1/2 to 4P3/2
  data["G0t3"] = np.array([0.379, 0.392, 0.394, 0.392, 0.392, 0.394, 0.389,
                           0.370, 0.326, 0.275, 0.221])
  # 2P1/2 to 4P5/2
  data["G0t4"] = np.array([0.242, 0.242, 0.242, 0.242, 0.247, 0.255, 0.253,
                           0.241, 0.210, 0.175, 0.141])
  # 2P3/2 to 4P3/2
  data["G1t3"] = np.array([0.401, 0.462, 0.479, 0.489, 0.501, 0.513, 0.509,
                           0.486, 0.428, 0.359, 0.289])
  # 2P3/2 to 4P5/2
  data["G1t4"] = np.array([1.138, 1.136, 1.131, 1.117, 1.111, 1.114, 1.098,
                           1.043, 0.919, 0.774, 0.629])
  # 4P1/2 to 4P3/2
  data["G2t3"] = np.array([0.574, 0.627, 0.643, 0.672, 0.792, 0.984, 1.073,
                           1.128, 1.118, 1.048, 0.919])
  # 4P1/2 to 4P5/2
  data["G2t4"] = np.array([0.660, 0.672, 0.685, 0.722, 0.836, 0.969, 1.005,
                           1.007, 0.975, 0.944, 0.911])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    imin = 1
    imax = 9
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 1000000,
                        p0 = (1, 100., 1., 1., 1., 0., 1.),
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
    code += get_code("CII", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.xlim(0., 1.e5)
    pl.savefig("tmp/CII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
