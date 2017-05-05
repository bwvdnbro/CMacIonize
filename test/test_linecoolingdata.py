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
# @file test_linecoolingdata.py
#
# @brief Unit test for the Python LineCoolingData module.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# Load the module.
import load_module as pycmi
# Load the sys module (for sys.exit())
import sys

import numpy as np

##
# @brief Unit test for the Python LineCoolingData module.
##
def main():
  file = open("linestr_testdata.txt")

  linecoolingdata = pycmi.liblinecoolingdata.LineCoolingData()

  lines = file.readlines()
  for line in lines:
    data = line.split()
    T = float(data[0])
    # convert to m^-3
    ne = float(data[1])*1.e6
    abundances = [float(ab) for ab in data[2:14]]
    fresults = [float(res) for res in data[14:]]

    results = linecoolingdata.linestr(np.array([T]), np.array([ne]), abundances)

    if len(results) != len(fresults):
      print "Error: expected {a} values, got {b}...".format(a = len(fresults),
                                                            b = len(results))
      sys.exit(1)
    names = ["c6300", "c9405", "c6312", "c33mu", "c19mu", "c3729", "c3727",
             "c7330", "c4363", "c5007", "c52mu", "c88mu", "c5755", "c6584",
             "c4072", "c6717", "c6725", "c3869", "cniii57", "cneii12",
             "cneiii15", "cnii122", "cii2325", "ciii1908", "coii7325",
             "csiv10"]
    for i in range(len(names)):
      a = results[names[i]][0]
      # convert from erg s^-1 to J s^-1
      b = fresults[i]*1.e-7
      if abs(a-b) > 1.e-14*abs(a+b):
        print "Error: {a} != {b} ({rel})!".format(a = a, b = b,
                                                  rel = abs(a-b)/abs(a+b))
        sys.exit(1)

  sys.exit(0)

# make sure the main function is executed
if __name__ == "__main__":
  main()
