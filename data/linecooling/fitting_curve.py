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
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##


##
# @brief Fitting curve.
#
# @param T Temperature value (in K).
# @param norm Value at 10,000 K, used to normalize the curve.
# @param exponent Value of the exponent.
# @return Value of the fitting curve.
##
def fitting_curve(T, norm, exponent):
  T4 = T * 1.e-4
  return norm * T4**exponent

##
# @brief Print the fit variables to the stdout.
#
# @param norm Value at 10,000 K, used to normalize the curve.
# @param exponent Value of the exponent.
##
def print_fit_variables(norm, exponent):
  print "Norm:", norm
  print "Exponent:", exponent

##
# @brief Initialize an empty dictionary to store the fitting variables in.
#
# @return Empty dictionary to store the fitting variables.
##
def initialize_data_values():
  return {"om": "", "ome": ""}

##
# @brief Append the given fit variables to the given dictionary.
#
# @param data_values Dictionary to append to.
# @param norm Value at 10,000 K, used to normalize the curve.
# @param exponent Value of the exponent.
##
def append_data_values(data_values, norm, exponent):
  data_values["om"] += "{value},".format(value = norm)
  data_values["ome"] += "{value},".format(value = exponent)

##
# @brief Print the given dictionary to the stdout.
#
# @param data_values Dictionary with fit variables.
##
def print_data_values(data_values):
  print "values omega:"
  print data_values["om"]
  print "values omega exponent:"
  print data_values["ome"]

##
# @brief Get the C++ code string for the given element and transition,
# containing the given fit variables.
#
# @param element IonName.
# @param transition TransitionName.
# @param norm Value at 10,000 K, used to normalize the curve.
# @param exponent Value of the exponent.
# @return C++ code string.
##
def get_code(element, transition, norm, exponent):
  code = "_collision_strength[{element}][{transition}] = {value};\n".format(
           element = element, transition = transition, value = norm)
  code += \
      "_collision_strength_exponent[CIII][{transition}] = {value};\n".format(
        element = element, transition = transition, value = exponent)
  return code
