
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
# @file fitting_curve.py
#
# @brief General fitting curve used to fit all data.
#
# The general shape of the curve was inspired by Burgess, A. & Tully, J. A.
# 1992, AAS, 254, 436 (http://adsabs.harvard.edu/abs/1992A&A...254..436B).
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np

##
# @brief Fitting curve.
#
# The fitting curve has the form
# \f[
#   AT^B \left( C + \frac{D}{T} + ET + F \log(T) + GT^H \right).
# \f]
#
# @param T Temperature value (in K).
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @param H Parameter \f$H\f$.
# @return Value of the fitting curve.
##
def fitting_curve(T, A, B, C, D, E, F, G, H):
  return A * T**B * (C + D / T + E * T + F * np.log(T) + G * T**H)

##
# @brief Print the fit variables to the stdout.
#
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @param H Parameter \f$H\f$.
##
def print_fit_variables(A, B, C, D, E, F, G, H):
  print "A, B, C, D, E, F, G, H:", A, B, C, D, E, F, G, H

##
# @brief Initialize an empty dictionary to store the fitting variables in.
#
# @return Empty dictionary to store the fitting variables.
##
def initialize_data_values():
  return {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": "",
          "H": ""}

##
# @brief Append the given fit variables to the given dictionary.
#
# @param data_values Dictionary to append to.
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @param H Parameter \f$H\f$.
##
def append_data_values(data_values, A, B, C, D, E, F, G, H):
  data_values["A"] += "{value},".format(value = A)
  data_values["B"] += "{value},".format(value = B)
  data_values["C"] += "{value},".format(value = C)
  data_values["D"] += "{value},".format(value = D)
  data_values["E"] += "{value},".format(value = E)
  data_values["F"] += "{value},".format(value = F)
  data_values["G"] += "{value},".format(value = G)
  data_values["H"] += "{value},".format(value = H)

##
# @brief Print the given dictionary to the stdout.
#
# @param data_values Dictionary with fit variables.
##
def print_data_values(data_values):
  print "values A:"
  print data_values["A"]
  print "values B:"
  print data_values["B"]
  print "values C:"
  print data_values["C"]
  print "values D:"
  print data_values["D"]
  print "values E:"
  print data_values["E"]
  print "values F:"
  print data_values["F"]
  print "values G:"
  print data_values["G"]
  print "values H:"
  print data_values["H"]

##
# @brief Get the C++ code string for the given element and transition,
# containing the given fit variables.
#
# @param element IonName.
# @param transition TransitionName.
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @param H Parameter \f$H\f$.
# @return C++ code string.
##
def get_code(element, transition, A, B, C, D, E, F, G, H):
  code  = "_collision_strength[{element}][{transition}][0] = {value};\n".format(
           element = element, transition = transition, value = A)
  code += "_collision_strength[{element}][{transition}][1] = {value};\n".format(
           element = element, transition = transition, value = B)
  code += "_collision_strength[{element}][{transition}][2] = {value};\n".format(
           element = element, transition = transition, value = C)
  code += "_collision_strength[{element}][{transition}][3] = {value};\n".format(
           element = element, transition = transition, value = D)
  code += "_collision_strength[{element}][{transition}][4] = {value};\n".format(
           element = element, transition = transition, value = E)
  code += "_collision_strength[{element}][{transition}][5] = {value};\n".format(
           element = element, transition = transition, value = F)
  code += "_collision_strength[{element}][{transition}][6] = {value};\n".format(
           element = element, transition = transition, value = G)
  code += "_collision_strength[{element}][{transition}][7] = {value};\n".format(
           element = element, transition = transition, value = H)
  code += "_collision_strength[{element}][{transition}][0] = {value};\n".format(
           element = element, transition = transition, value = A)
  return code
