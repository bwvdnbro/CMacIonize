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
# @file gamma_SII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Tayal, S. S. & Zatsarinny, O. 2010, ApJS, 188, 32
# (http://adsabs.harvard.edu/abs/2010ApJS..188...32T).
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

  # data from Tayal & Zatsarinny (2010), table 5
  T = np.array([5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000.,
                50000., 70000., 100000.])
  # 4S3/2 to 2D3/2
  data["G0t1"] = np.array([2.66e+0, 2.62e+0, 2.56e+0, 2.48e+0, 2.41e+0,
                           2.35e+0, 2.30e+0, 2.20e+0, 2.10e+0, 1.93e+0,
                           1.71e+0])
  # 4S3/2 to 2D5/2
  data["G0t2"] = np.array([3.98e+0, 3.91e+0, 3.83e+0, 3.71e+0, 3.61e+0,
                           3.52e+0, 3.44e+0, 3.30e+0, 3.16e+0, 2.90e+0,
                           2.58e+0])
  # 2D3/2 to 2D5/2
  data["G1t2"] = np.array([7.32e+0, 7.13e+0, 6.89e+0, 6.58e+0, 6.35e+0,
                           6.17e+0, 6.02e+0, 5.76e+0, 5.51e+0, 5.06e+0,
                           4.48e+0])
  # 2P1/2 to 2P3/2
  data["G3t4"] = np.array([1.76e+0, 1.78e+0, 1.80e+0, 1.82e+0, 1.84e+0,
                           1.85e+0, 1.86e+0, 1.84e+0, 1.80e+0, 1.70e+0,
                           1.54e+0])
  # 4S3/2 to 2P1/2
  data["G0t3"] = np.array([6.86e-1, 6.94e-1, 7.04e-1, 7.17e-1, 7.27e-1,
                           7.33e-1, 7.36e-1, 7.30e-1, 7.15e-1, 6.70e-1,
                           5.97e-1])
  # 4S3/2 to 2P3/2
  data["G0t4"] = np.array([1.38e+0, 1.39e+0, 1.42e+0, 1.44e+0, 1.46e+0,
                           1.47e+0, 1.48e+0, 1.47e+0, 1.43e+0, 1.34e+0,
                           1.19e+0])
  # 2D3/2 to 2P1/2
  data["G1t3"] = np.array([1.48e+0, 1.48e+0, 1.47e+0, 1.48e+0, 1.48e+0,
                           1.49e+0, 1.50e+0, 1.50e+0, 1.50e+0, 1.49e+0,
                           1.45e+0])
  # 2D3/2 to 2P3/2
  data["G1t4"] = np.array([2.40e+0, 2.40e+0, 2.39e+0, 2.39e+0, 2.38e+0,
                           2.38e+0, 2.38e+0, 2.36e+0, 2.33e+0, 2.25e+0,
                           2.13e+0])
  # 2D5/2 to 2P1/2
  data["G2t3"] = np.array([1.79e+0, 1.78e+0, 1.78e+0, 1.77e+0, 1.76e+0,
                           1.76e+0, 1.75e+0, 1.73e+0, 1.70e+0, 1.63e+0,
                           1.53e+0])
  # 2D5/2 to 2P3/2
  data["G2t4"] = np.array([4.07e+0, 4.07e+0, 4.06e+0, 4.06e+0, 4.07e+0,
                           4.08e+0, 4.08e+0, 4.08e+0, 4.06e+0, 3.97e+0,
                           3.83e+0])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  for key in sorted(data):
    # we force the curve to go through the value at 10,000 K using a global
    # variable
    norm = data[key][2]
    # we start by fitting to the full data set
    imin = 0
    imax = len(T)
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[key][imin:imax])
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum( (data[key][imin:imax] - fitting_curve(T[imin:imax], *A))**2 )
    # if xi2 is too large: shrink the fitting interval and try again
    while xi2 > 1.e-5:
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
    pl.savefig("tmp/SII_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print_data_values(data_values)
