#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file test_densitygrid.py
#
# @brief Unit test for the Python DensityGrid module.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# Load the module.
import load_module as pycmi
# Load the sys module (for sys.exit())
import sys

import numpy as np
import h5py

##
# @brief Unit test for the Python DensityGrid module.
##
def main():
  densitygrid = pycmi.libdensitygrid.DensityGrid("python_test.hdf5")

  file = h5py.File("python_test.hdf5", 'r')
  fcoords = np.array(file["/PartType0/Coordinates"])
  fbox = np.array(file["/Header"].attrs["BoxSize"])

  ncell = len(fcoords)

  if not densitygrid.get_number_of_cells() == ncell:
    print "Error: wrong number of cells ({a}, expected: {b})!".format(
      a = densitygrid.get_number_of_cells(), b = ncell)
    sys.exit(1)

  box = densitygrid.get_box()
  for i in range(3):
    if not box["sides"][i] == fbox[i]:
      print "Error: wrong box size ([{a}, {b}, {c}], " \
            "expected [{d}, {e}, {f}])!".format(
        a = box["sides"][0], b = box["sides"][1], c = box["sides"][2],
        d = fbox[0], e = fbox[1], f = fbox[2])
      sys.exit(1)
  if not box["units"] == "m":
    print "Error: wrong unit (\"{a}\", expected \"{b}\")!".format(
      a = box["units"], b = "m")
    sys.exit(1)

  coords = densitygrid.get_coordinates()
  if not len(coords["values"]) == ncell:
    print "Error: wrong number of coordinates ({a}, expected {b})!".format(
      a = len(coords["values"]), b = ncell)
    sys.exit(1)
  if not len(coords["values"][0]) == 3:
    print "Error: wrong shape (({a}, {b}), expected ({c}, {d}))!".format(
      a = ncell, b = len(coords["values"][0]), c = ncell, d = 3)
    sys.exit(1)
  for i in range(ncell):
    for j in range(3):
      if not (coords["values"][i][j] - box["origin"][j]) == fcoords[i][j]:
        print "Error: wrong coordinates ([{a}, {b}, {c}], " \
              "expected [{d}, {e}, {f}])!".format(
          a = coords["values"][i][0] - box["origin"][0],
          b = coords["values"][i][1] - box["origin"][1],
          c = coords["values"][i][2] - box["origin"][2],
          d = fcoords[i][0], e = fcoords[i][1], f = fcoords[i][2])
        sys.exit(1)
  if not coords["units"] == "m":
    print "Error: wrong unit (\"{a}\", expected \"{b}\")!".format(
      a = coords["units"], b = "m")
    sys.exit(1)

  units = {"NumberDensity": "m^-3", "Temperature": "K", "NeutralFractionH": "",
           "NeutralFractionHe": "", "NeutralFractionC+": "",
           "NeutralFractionC++": "", "NeutralFractionN": "",
           "NeutralFractionN+": "", "NeutralFractionN++": "",
           "NeutralFractionO": "", "NeutralFractionO+": "",
           "NeutralFractionNe": "", "NeutralFractionNe+": "",
           "NeutralFractionS+": "", "NeutralFractionS++": "",
           "NeutralFractionS+++": ""}
  for variable in units:
    data = densitygrid.get_variable(variable)
    fdata = np.array(file["/PartType0/{variable}".format(variable = variable)])
    if not len(data["values"]) == len(fdata):
      print "Error: wrong number of values ({a}, expected {b})!".format(
        a = len(data["values"]), b = len(fdata))
      sys.exit(1)
    for i in range(len(fdata)):
      if not data["values"][i] == fdata[i]:
        print "Error: wrong {variable} value ({a}, expected {b})!".format(
          variable = variable, a = data["values"][i], b = fdata[i])
        sys.exit(1)
    if not data["units"] == units[variable]:
      print "Error: wrong unit (\"{a}\", expected \"{b}\")!".format(
        a = data["units"], b = units[variable])
      sys.exit(1)

  sys.exit(0)

# make sure the main function is executed
if __name__ == "__main__":
  main()
