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
# @file gamma_OI.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Zatsarinny, O. & Tayal, S. S. 2003, ApJS, 148, 575
# (http://adsabs.harvard.edu/abs/2003ApJS..148..575Z) (ZT) and Berrington, K. A.
# 1988, JPhB, 21, 1083 (http://adsabs.harvard.edu/abs/1988JPhB...21.1083B) (B).
# We also use a formula from Lennon, D. J. & Burke, V. M. 1994, A&AS, 103, 273
# (http://adsabs.harvard.edu/abs/1994A%26AS..103..273L).
#
# Note that Bell, K. L., Berrington, K. A. & Thomas, M. R. J. 1998, MNRAS, 293,
# L83 (http://adsabs.harvard.edu/abs/1998MNRAS.293L..83B) provide more accurate
# low temperature fine structure collision strengths, but these do not go up to
# 10,000 K. Values above 10,000 K do not exist in literature.
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

# dictionary that links abbreviated transition names to the full names used in
# LineCoolingData
transitions = {
  "G0t1": ["TRANSITION_0_to_1", "B"],
  "G0t2": ["TRANSITION_0_to_2", "B"],
  "G0t3": ["TRANSITION_0_to_3", "ZT"],
  "G0t4": ["TRANSITION_0_to_4", "ZT"],
  "G1t2": ["TRANSITION_1_to_2", "B"],
  "G1t3": ["TRANSITION_1_to_3", "ZT"],
  "G1t4": ["TRANSITION_1_to_4", "ZT"],
  "G2t3": ["TRANSITION_2_to_3", "ZT"],
  "G2t4": ["TRANSITION_2_to_4", "ZT"],
  "G3t4": ["TRANSITION_3_to_4", "ZT"]
}

# index of the 10,000 K element in the collision strength array
tKindex = {"ZT": 9, "B": 7}

##
# @brief Fitting curve.
#
# @param T Temperature value (in K).
# @param A Value of the exponent.
##
def fitting_curve(T, A):
  T4 = T * 1.e-4
  return norm * T4**A

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
  T = {}
  data = {}

  # data from Berrington (1988), table 2
  T["B"] = np.array([50., 100., 200., 500., 1000., 2000., 5000., 10000.])
  # 3P2 to 3P1
  data["G0t1"] = np.array([1.08e-3, 1.74e-3, 2.82e-3, 5.69e-3, 9.72e-3,
                           1.7e-2, 4.52e-2, 0.106])
  # 3P2 to 3P0
  data["G0t2"] = np.array([5.42e-4, 8.22e-4, 1.28e-3, 2.51e-3, 4.11e-3,
                           6.78e-3, 1.53e-2, 3.21e-2])
  # 3P1 to 3P0
  data["G1t2"] = np.array([4.05e-5, 6.62e-5, 1.16e-4, 2.86e-4, 6.36e-4,
                           1.67e-3, 9.12e-3, 2.83e-2])

  # data from Zatsarinny & Tayal (2003)
  T["ZT"] = np.array([1000., 1500., 2000., 2500., 3000., 4000., 5000., 6000.,
                      8000., 10000., 12000., 14000., 16000., 20000., 25000.,
                      30000., 35000., 40000., 50000., 60000.])
  # 3P? to 1D2
  GPtD = np.array([3.17e-2, 4.49e-2, 5.87e-2, 7.29e-2, 8.74e-2, 1.17e-1,
                   1.46e-1, 1.76e-1, 2.35e-1, 2.93e-1, 3.48e-1, 4.01e-1,
                   4.50e-1, 5.41e-1, 6.41e-1, 7.26e-1, 8.00e-1, 8.65e-1,
                   9.75e-1, 1.06])
  # 3P? to 1S0
  GPtS = np.array([3.14e-3, 3.60e-3, 4.49e-3, 5.73e-3, 7.21e-3, 1.06e-2,
                   1.43e-2, 1.80e-2, 2.53e-2, 3.23e-2, 3.88e-2, 4.49e-2,
                   5.05e-2, 6.07e-2, 7.17e-2, 8.11e-2, 8.92e-2, 9.62e-2,
                   1.08e-1, 1.17e-1])
  # 1D2 to 1S0
  data["G3t4"] = np.array([1.01e-2, 1.74e-2, 2.44e-2, 3.09e-2, 3.67e-2,
                           4.69e-2, 5.56e-2, 6.33e-2, 7.67e-2, 8.83e-2,
                           9.86e-2, 1.08e-1, 1.16e-1, 1.32e-1, 1.50e-1,
                           1.66e-1, 1.81e-1, 1.95e-1, 2.21e-1, 2.45e-1])

  # formula (1) in Lennon & Burke (1994)
  data["G0t3"] = 5.*GPtD/9.
  data["G0t4"] = 5.*GPtS/9.
  data["G1t3"] = GPtD/3.
  data["G1t4"] = GPtS/3.
  data["G2t3"] = GPtD/9.
  data["G2t4"] = GPtS/9.

  # initialize the strings for code and value output
  code = ""
  values_om = ""
  values_ome = ""
  # do the curve fitting
  for key in sorted(data):
    dkey = transitions[key][1]
    # we force the curve to go through the value at 10,000 K using a global
    # variable
    norm = data[key][tKindex[dkey]]
    # we start by fitting to the full data set
    imin = 0
    imax = len(T[dkey])
    # fit the curve
    A,_ = opt.curve_fit(fitting_curve, T[dkey][imin:imax], data[key][imin:imax])
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = \
      sum( (data[key][imin:imax] - fitting_curve(T[dkey][imin:imax], *A))**2 )
    # if xi2 is too large: shrink the fitting interval and try again
    while xi2 > 1.e100:
      imin += 1
      imax -= 1
      # if the interval becomes too small, we bail out
      if imax - imin < 2:
        break
      A,_ = opt.curve_fit(fitting_curve, T[dkey][imin:imax],
                          data[key][imin:imax])
      xi2 = \
        sum( (data[key][imin:imax] - fitting_curve(T[dkey][imin:imax], *A))**2 )
    # output some info
    print "Transition:", key
    print "Gamma:", norm
    print "exponent:", A
    print "convergence:", xi2
    print "validity: [", T[dkey][imin], ",", T[dkey][imax-1], "]"
    # write the fitting code for this transition
    code += "_collision_strength[OI][{transition}] = {value};\n".format(
      transition = transitions[key][0], value = norm)
    code += \
      "_collision_strength_exponent[OI][{transition}] = {value};\n".format(
        transition = transitions[key][0], value = A[0])
    # add the values to the list strings
    values_om += "{value},".format(value = norm)
    values_ome += "{value},".format(value = A[0])

    # plot the data and fit for visual comparison
    Trange = np.logspace(3., 5., 100)
    pl.plot(T[dkey], data[key], "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.savefig("tmp/OI_{key}.png".format(key = key))
    pl.close()

  # output the code to put into the LineCoolingData constructor
  print "code:"
  print code
  # output the values to put in atom4.dat in Kenny's code (to update the
  # reference values for the unit tests)
  print "values omega:"
  print values_om
  print "values omega exponent:"
  print values_ome
