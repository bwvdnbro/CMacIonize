#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file write_LambertWtest.py
#
# @brief Write a reference file that can be used to test our custom
# implementations of the LambertW 0 and -1 branches.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import scipy.special.lambertw as lambertw

# we add a very small extra offset to the first point, because otherwise we
# would end up with the wrong value for LambertW_-1 in that point...
x = np.arange(-1.+1.e-15, 0., 0.001) / np.exp(1.)
f0 = lambertw(x, 0).real
fm1 = lambertw(x, -1).real

file = open("test_lambertw_input.txt", "w")
file.write("# file written by write_LambertWtest.py\n")
file.write("# x\tLambertW_0\tLambertW_-1\n")
for i in range(len(x)):
  file.write("{0}\t{1}\t{2}\n".format(x[i], f0[i], fm1[i]))
file.close()
