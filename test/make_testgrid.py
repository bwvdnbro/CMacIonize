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
# @file make_testgrid.py
#
# @brief Write the testgrid.txt ASCII test file to test the
# AsciiFileDensityFunction class.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

file = open("testgrid.txt", 'w')

for i in range(8):
  for j in range(8):
    for k in range(8):
      x = (i+0.5)*0.125
      y = (j+0.5)*0.125
      z = (k+0.5)*0.125
      file.write("{x}\t{y}\t{z}\t{rho}\n".format(x = x, y = y, z = z, rho = 1.))
