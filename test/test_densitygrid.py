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

##
# @brief Unit test for the Python DensityGrid module.
##
def main():
  densitygrid = pycmi.libdensitygrid.DensityGrid("python_test.hdf5")

  if not densitygrid.get_number_of_cells() == 64:
    print "Error: number of cells {a} (expected: {b})!".format(
      a = densitygrid.get_number_of_cells(), b = 64)
    sys.exit(1)

  sys.exit(0)

# make sure the main function is executed
if __name__ == "__main__":
  main()
