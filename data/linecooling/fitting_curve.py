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
#   T^{A+1} \left( B + \frac{C}{T} + D \log(T) +
#                  E T \left( 1 + (F - 1) T^G \right) \right).
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
  return T**(A + 1.) * \
         (B + C / T + D * np.log(T) + E * T * (1. + (F - 1.) * T**G))

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
  J[:, 1] = T**(A + 1.)
  J[:, 2] = T**A
  J[:, 3] = T**(A + 1.) * np.log(T)
  J[:, 4] = T**(A + 2.) + (F - 1.) * T**(A + G + 2.)
  J[:, 5] = E * T**(A + G + 2.)
  J[:, 6] = np.log(T) * E * (F - 1.) * T**(A + G + 2.)
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
  return {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": "",
          "table": []}

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
  data_values["A"] += "{value:.3e},".format(value = A)
  data_values["B"] += "{value:.3e},".format(value = B)
  data_values["C"] += "{value:.3e},".format(value = C)
  data_values["D"] += "{value:.3e},".format(value = D)
  data_values["E"] += "{value:.3e},".format(value = E)
  data_values["F"] += "{value:.3e},".format(value = F)
  data_values["G"] += "{value:.3e},".format(value = G)
  data_values["table"].append({})
  data_values["table"][-1]["A"] = "{0:.3e}".format(A)
  data_values["table"][-1]["B"] = "{0:.3e}".format(B)
  data_values["table"][-1]["C"] = "{0:.3e}".format(C)
  data_values["table"][-1]["D"] = "{0:.3e}".format(D)
  data_values["table"][-1]["E"] = "{0:.3e}".format(E)
  data_values["table"][-1]["F"] = "{0:.3e}".format(F)
  data_values["table"][-1]["G"] = "{0:.3e}".format(G)

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
  print "table:"
  print r"""$\begin{matrix}
a \\
b \\
c \\
d \\
e \\
f \\
g \\
\end{matrix}$"""
  for i in range(5):
    print "&\n\\begin{{tabular}}{{l}}{a} \\\\{b} \\\\{c} \\\\{d} \\\\".format(
      a = data_values["table"][i]["A"], b = data_values["table"][i]["B"],
      c = data_values["table"][i]["C"], d = data_values["table"][i]["D"])
    print "{e} \\\\{f} \\\\{g} \\\\\\end{{tabular}}".format(
      e = data_values["table"][i]["E"], f = data_values["table"][i]["F"],
      g = data_values["table"][i]["G"])
  print r"""\vspace{12pt}
\\
&
$\begin{matrix}
a \\
b \\
c \\
d \\
e \\
f \\
g \\
\end{matrix}$"""
  for i in range(5, 10):
    print "&\n\\begin{{tabular}}{{l}}{a} \\\\{b} \\\\{c} \\\\{d} \\\\".format(
      a = data_values["table"][i]["A"], b = data_values["table"][i]["B"],
      c = data_values["table"][i]["C"], d = data_values["table"][i]["D"])
    print "{e} \\\\{f} \\\\{g} \\\\\\end{{tabular}}".format(
      e = data_values["table"][i]["E"], f = data_values["table"][i]["F"],
      g = data_values["table"][i]["G"])

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
  code  = "_collision_strength[{e}][{t}][0] = {v:.3e};\n".format(
           e = element, t = transition, v = A)
  code += "_collision_strength[{e}][{t}][1] = {v:.3e};\n".format(
           e = element, t = transition, v = B)
  code += "_collision_strength[{e}][{t}][2] = {v:.3e};\n".format(
           e = element, t = transition, v = C)
  code += "_collision_strength[{e}][{t}][3] = {v:.3e};\n".format(
           e = element, t = transition, v = D)
  code += "_collision_strength[{e}][{t}][4] = {v:.3e};\n".format(
           e = element, t = transition, v = E)
  code += "_collision_strength[{e}][{t}][5] = {v:.3e};\n".format(
           e = element, t = transition, v = F)
  code += "_collision_strength[{e}][{t}][6] = {v:.3e};\n".format(
           e = element, t = transition, v = G)
  return code

##
# @brief Round the parameter values to 4 significant digits.
#
# @param A Parameter \f$A\f$.
# @param B Parameter \f$B\f$.
# @param C Parameter \f$C\f$.
# @param D Parameter \f$D\f$.
# @param E Parameter \f$E\f$.
# @param F Parameter \f$F\f$.
# @param G Parameter \f$G\f$.
# @return Rounded parameter values.
##
def round_parameters(A, B, C, D, E, F, G):
  A = float("{0:.3e}".format(A))
  B = float("{0:.3e}".format(B))
  C = float("{0:.3e}".format(C))
  D = float("{0:.3e}".format(D))
  E = float("{0:.3e}".format(E))
  F = float("{0:.3e}".format(F))
  G = float("{0:.3e}".format(G))
  return A, B, C, D, E, F, G
