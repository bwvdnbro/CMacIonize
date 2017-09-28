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
# @file gamma_NIII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Blum, R. D. & Pradhan, A. K. 1992, ApJS, 80, 425
# (http://adsabs.harvard.edu/abs/1992ApJS...80..425B).
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

# main function: computes fits to the data of Blum & Pradhan (1992) and plots
# the data and fits for visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  # data from Blum & Pradhan (1992), table 3, 1 to 2 transition
  T = np.array([1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000.,
                10000., 11000., 12000., 13000., 14000., 15000., 16000., 18000.,
                20000., 22000., 24000., 26000., 28000., 30000., 32000., 34000.,
                36000., 38000., 40000.])
  data = np.array([1.2896, 1.2857, 1.2871, 1.2976, 1.3163, 1.3403, 1.3667,
                   1.3937, 1.42, 1.4454, 1.4694, 1.4922, 1.5139, 1.5345, 1.5541,
                   1.5729, 1.6086, 1.6421, 1.674, 1.7045, 1.7339, 1.7622,
                   1.7895, 1.8156, 1.8407, 1.8646, 1.8874, 1.909])

  # initialize the strings for code and value output
  code = ""
  data_values = initialize_data_values()
  # do the curve fitting
  imin = 0
  imax = len(T)
  # fit the curve
  A,_ = opt.curve_fit(fitting_curve, T[imin:imax], data[imin:imax],
                      maxfev = 100000,
                      p0 = (0., 1., 10., 1., 0.01, 0.01, 0.),
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
  code += get_code("NIII", "REMOVE_THIS_BRACKET", *A)
  # add the values to the list strings
  append_data_values(data_values, *A)

  # plot the data and fit for visual comparison
  Trange = np.logspace(3., 5., 100)
  pl.plot(T, data, "k.")
  pl.plot(Trange, fitting_curve(Trange, *A), "r-")
  pl.savefig("tmp/NIII_{key}.png".format(key = "G0t1"))
  pl.close()

  # save the plot values in separate files
  dfile = open("tmp/NIII_{key}_data.txt".format(key = "G0t1"), "w")
  for i in range(len(T)):
    dfile.write("{T}\t{data}\n".format(T = T[i], data = data[i]))
  dfile.close()
  ffile = open("tmp/NIII_{key}_fit.txt".format(key = "G0t1"), "w")
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
