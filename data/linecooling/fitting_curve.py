
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
#   T^A \left( B + \frac{C}{T} + DT + E \log(T) + FT^G \right).
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
# @return Value of the fitting curve.
##
def fitting_curve(T, A, B, C, D, E, F, G):
  return T**A * (B + C / T + D * T + E * np.log(T) + F * T**G)

##
# @brief Jacobian of the fitting curve.
#
# @param T Temperature value (in K).
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @return Value of the Jacobian of the fitting curve.
##
def jacobian_fitting_curve(T, A, B, C, D, E, F, G):
  J = np.zeros((len(T), 7))
  J[:, 0] = np.log(T) * fitting_curve(T, A, B, C, D, E, F, G)
  J[:, 1] = T**A
  J[:, 2] = T**(A - 1.)
  J[:, 3] = T**(A + 1.)
  J[:, 4] = T**A * np.log(T)
  J[:, 5] = T**(A + G)
  J[:, 6] = np.log(T) * F * T**(A + G)
  return J

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
##
def print_fit_variables(A, B, C, D, E, F, G):
  print "A, B, C, D, E, F, G:", A, B, C, D, E, F, G

##
# @brief Initialize an empty dictionary to store the fitting variables in.
#
# @return Empty dictionary to store the fitting variables.
##
def initialize_data_values():
  return {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": ""}

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
##
def append_data_values(data_values, A, B, C, D, E, F, G):
  data_values["A"] += "{value},".format(value = A)
  data_values["B"] += "{value},".format(value = B)
  data_values["C"] += "{value},".format(value = C)
  data_values["D"] += "{value},".format(value = D)
  data_values["E"] += "{value},".format(value = E)
  data_values["F"] += "{value},".format(value = F)
  data_values["G"] += "{value},".format(value = G)

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
# @return C++ code string.
##
def get_code(element, transition, A, B, C, D, E, F, G):
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
  return code
