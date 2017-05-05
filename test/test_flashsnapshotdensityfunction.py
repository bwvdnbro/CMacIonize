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
# @file test_flashsnapshotdensityfunction.py
#
# @brief Unit test for the Python FLASHSnapshotDensityFunction module.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# Load the module.
import load_module as pycmi
# alias the libflashsnapshotdensityfunction for more readable code
flashlib = pycmi.libflashsnapshotdensityfunction
# Load the sys module (for sys.exit())
import sys

##
# @brief Unit test for the Python FLASHSnapshotDensityFunction module.
##
def main():
  flashfunc = flashlib.FLASHSnapshotDensityFunction("FLASHtest.hdf5")

  temp = flashfunc.get_temperature(0., 0., 0.)
  if temp != 4000.:
    print "Error: temperature incorrect ({val}, expected {expval})!".format(
      val = temp, expval = 4000.)
    sys.exit(1)

  sys.exit(0)

# make sure the main function is executed
if __name__ == "__main__":
  main()
