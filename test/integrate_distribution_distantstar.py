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
# @file integrate_distribution_distantstar.py
#
# @brief Script to calculate the expected value of the 2D distribution function
# for the positions generated in testDistantStarContinuousPhotonSource.cpp.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# we really need Scientific Python's power here
import scipy.integrate as quad

##
# @brief Distribution function describing the distribution of positions where
# photons emitted by an isotropic source on a position (2., 2., 2.) hit the part
# of the plane z = 0 between x = [0., 1.] and y = [0., 1.].
#
# @param x X-position in the range [0., 1.]
# @param y Y-position in the range [0., 1.]
# @return Value of the non-normalized distribution function.
##
def distribution(x, y):
  return 1./((x-2.)**2 + (y-2.)**2 + 1.)**1.5

##
# @brief Integrand of the expectation value integral.
#
# This is just x * distribution. We limit ourselves to a single coordinate; the
# other coordinates have the same expectation value due to symmetry.
#
# @param x X-position in the range [0., 1.]
# @param y Y-position in the range [0., 1.]
# @return Value of the expectation value integrand.
##
def avg_int(x, y):
  return x * distribution(x, y)

##
# @brief Print the expectation value of the distribution function to the
# standard output.
##
def main():
  # The distribution function was not normalized, so we first calculate the
  # normalization constant (which is 1./norm)
  # scipy.integrate.dblquad requires two functions as x limits. Since our x
  # limits do not depend on y, we provide trivial lambda function versions.
  norm = quad.dblquad(distribution, 0., 1., lambda x: 0., lambda x: 1.)[0]
  # Now we calculate the expectation value for a single coordinate.
  avg = quad.dblquad(avg_int, 0., 1., lambda x: 0., lambda x: 1.)[0]
  # And we print the result.
  print "Expected value:", avg/norm

# Make sure main() is called.
if __name__ == "__main__":
  main()
