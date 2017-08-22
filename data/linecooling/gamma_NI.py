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
# @file gamma_NI.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Tayal, S. S. 2000, ADNDT, 76, 191
# (http://adsabs.harvard.edu/abs/2000ADNDT..76..191T).
#
# Note that for consistency we assume the same level order as in Froese Fischer,
# C. & Tachiev, G. 2004, ADNDT, 87, 1
# (http://adsabs.harvard.edu/abs/2004ADNDT..87....1F), which means that some of
# the transitions do not correspond to the levels in Tayal (2000).
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

  # data from Tayal (2000), table III
  T = np.array([1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500., 5000.,
                6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000.,
                25000., 30000., 35000., 40000., 50000., 60000., 80000., 100000.,
                200000., 300000., 400000., 500000., 600000.])
  # 4S3/2 to 2D5/3
  data["G0t1"] = np.array([4.58e-4, 7.74e-4, 1.13e-3, 1.52e-3, 1.92e-3, 2.35e-3,
                           2.81e-3, 3.33e-3, 3.96e-3, 5.74e-3, 1.28e-2, 2.62e-2,
                           4.62e-2, 7.16e-2, 1.01e-1, 1.33e-1, 1.66e-1, 2.50e-1,
                           3.28e-1, 3.99e-1, 4.62e-1, 5.65e-1, 6.45e-1, 7.59e-1,
                           8.34e-1, 9.80e-1, 9.89e-1, 9.55e-1, 9.06e-1,
                           8.54e-1])
  # 4S3/2 to 2D3/2
  data["G0t2"] = np.array([3.05e-4, 5.16e-4, 7.54e-4, 1.01e-3, 1.28e-3, 1.57e-3,
                           1.87e-3, 2.22e-3, 2.64e-3, 3.83e-3, 8.54e-3, 1.75e-2,
                           3.08e-2, 4.77e-2, 6.73e-2, 8.86e-2, 1.11e-1, 1.66e-1,
                           2.19e-1, 2.66e-1, 3.08e-1, 3.77e-1, 4.30e-1, 5.06e-1,
                           5.56e-1, 6.54e-1, 6.59e-1, 6.36e-1, 6.04e-1,
                           5.69e-1])
  # 2D5/2 to 2D3/2
  data["G1t2"] = np.array([3.69, 3.92, 3.98, 3.95, 3.88, 3.78, 3.67, 3.56, 3.45,
                           3.27, 3.09, 3.24, 3.69, 4.37, 5.20, 6.11, 7.07, 9.45,
                           11.7, 13.7, 15.5, 18.4, 20.6, 23.4, 24.8, 24.1, 20.9,
                           18.1, 15.9, 14.2])
  # 2P1/2 to 2P3/2
  data["G3t4"] = np.array([1.60, 1.66, 1.68, 1.67, 1.65, 1.62, 1.58, 1.54, 1.51,
                           1.47, 1.57, 1.99, 2.68, 3.58, 4.58, 5.61, 6.63, 8.92,
                           10.8, 12.1, 13.1, 14.3, 14.8, 14.8, 14.3, 11.1, 8.86,
                           7.35, 6.29, 5.49])
  # 4S3/2 to 2P1/2
  data["G0t3"] = np.array([1.34e-3, 2.06e-3, 2.76e-3, 3.44e-3, 4.09e-3, 4.71e-3,
                           5.29e-3, 5.83e-3, 6.34e-3, 7.26e-3, 8.85e-3, 1.05e-2,
                           1.25e-2, 1.50e-2, 1.80e-2, 2.15e-2, 2.53e-2, 3.57e-2,
                           4.63e-2, 5.63e-2, 6.54e-2, 8.11e-2, 9.37e-2, 1.13e-1,
                           1.26e-1, 1.56e-1, 1.61e-1, 1.58e-1, 1.50e-1,
                           1.42e-1])
  # 4S3/2 to 2P3/2
  data["G0t4"] = np.array([2.68e-3, 4.11e-3, 5.51e-3, 6.87e-3, 8.18e-3, 9.41e-3,
                           1.06e-2, 1.17e-2, 1.27e-2, 1.45e-2, 1.77e-2, 2.10e-2,
                           2.50e-2, 3.00e-2, 3.61e-2, 4.30e-2, 5.07e-2, 7.14e-2,
                           9.25e-2, 1.13e-1, 1.31e-1, 1.62e-1, 1.87e-1, 2.25e-1,
                           2.51e-1, 3.12e-1, 3.23e-1, 3.15e-1, 3.01e-1,
                           2.85e-1])
  # 2D5/2 to 2P1/2
  data["G1t3"] = np.array([4.65e-2, 5.18e-2, 5.24e-2, 5.17e-2, 5.07e-2, 5.05e-2,
                           5.19e-2, 5.63e-2, 6.50e-2, 1.00e-1, 2.64e-1, 5.60e-1,
                           9.57e-1, 1.41, 1.89, 2.35, 2.80, 3.75, 4.48, 5.00,
                           5.38, 5.82, 6.00, 6.01, 5.83, 4.67, 3.86, 3.30, 2.89,
                           2.58])
  # 2D5/2 to 2P3/2
  data["G1t4"] = np.array([9.09e-2, 1.02e-1, 1.04e-1, 1.03e-1, 1.04e-1, 1.12e-1,
                           1.33e-1, 1.80e-1, 2.61e-1, 5.70e-1, 1.91, 4.15, 6.96,
                           10.0, 13.0, 15.9, 18.5, 23.8, 27.5, 30.0, 31.6, 33.0,
                           33.1, 31.8, 29.8, 21.7, 17.0, 14.1, 12.1, 10.6])
  # 2D3/2 to 2P1/2
  data["G2t3"] = np.array([2.99e-2, 3.34e-2, 3.42e-2, 3.44e-2, 3.54e-2, 3.97e-2,
                           5.13e-2, 7.53e-2, 1.18e-1, 2.77e-1, 9.67e-1, 2.11,
                           3.54, 5.08, 6.59, 8.00, 9.28, 11.8, 13.6, 14.7, 15.4,
                           16.0, 15.9, 15.2, 14.1, 10.1, 7.87, 6.49, 5.54,
                           4.85])
  # 2D3/2 to 2P3/2
  data["G2t4"] = np.array([6.17e-2, 6.88e-2, 6.98e-2, 6.89e-2, 6.79e-2, 6.85e-2,
                           7.25e-2, 8.24e-2, 1.01e-1, 1.74e-1, 5.02e-1, 1.07,
                           1.82, 2.66, 3.52, 4.36, 5.14, 6.80, 8.03, 8.91, 9.51,
                           10.2, 10.4, 10.3, 9.88, 7.69, 6.26, 5.30, 4.62,
                           4.10])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    imin = 7
    imax = 24
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax],
                        maxfev = 1000000,
                        p0 = (1., 1., 1., 1., 10., 1., 1.),
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
    code += get_code("NI", transitions[key], *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T, data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.xlim(0., 1.e5)
    pl.savefig("tmp/NI_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
