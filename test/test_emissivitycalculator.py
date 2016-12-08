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
# @file test_emissivitycalculator.py
#
# @brief Unit test for the Python EmissivityCalculator module.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# Load the module.
import load_module as pycmi
# Load the sys module (for sys.exit())
import sys

import numpy as np

##
# @brief Unit test for the Python EmissivityCalculator module.
##
def main():
  densitygrid = pycmi.libdensitygrid.DensityGrid("python_test.hdf5")

  abundances = {"helium": 0.1, "carbon": 2.2e-4, "nitrogen": 4.e-5,
                "oxygen": 3.3e-4, "neon": 5.e-5, "sulphur": 9.e-6}
  abundances = pycmi.libemissivitycalculator.Abundances(abundances)
  emissivitycalculator = \
      pycmi.libemissivitycalculator.EmissivityCalculator(abundances)

  emissivities = emissivitycalculator.get_emissivities(densitygrid)

  print emissivities

  sys.exit(0)

# make sure the main function is executed
if __name__ == "__main__":
  main()
