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
# @file emissivitycalculator.py
#
# @brief Python version of the EmissivityCalculator.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##
import liblinecoolingdata

# global object that is initialized once at the import of the module
lines = liblinecoolingdata.LineCoolingData()

##
# @brief Filter for cells for which we do not want to calculate emissivities.
#
# @param fractions Neutral and ionic fractions in the cell.
# @param temperature Temperature in the cell (in K).
# @return True if it is safe to calculate an emissivity for this cell.
##
def has_emissivity(fractions, temperature):
  return (fractions[0] < 0.2 and temperature > 3000.)

##
# @brief Get the emissivity values for a cell with the given ionic and neutral
# fractions, abundances, density, and temperature.
#
# @param fractions Neutral and ionic fractions in the cell: [H, He, [Cp1, Cp2],
# [N, Np1, Np2], [O, Op1], [Ne, Nep1], [Sp1, Sp2, Sp3]].
# @param abundances Abundances of all elements other than hydrogen (relative
# to hydrogen): [He, C, N, O, Ne, S].
# @param density Particle density in the cell (in m^-3).
# @param temperature Temperature in the cell (in K).
# @return Emissivity values in the cell (as a dictionary).
##
def get_emissivity(fractions, abundances, density, temperature):
  emissivities = {"Hbeta": 0.,
                  "HeI_5876": 0.,
                  "OI_6300": 0.}
  if(has_emissivity(fractions, temperature)):
    nhp = density * (1. - fractions[0])
    nhep = density * (1. - fractions[1]) * abundances[0]
    ne = nhp + nhep

    abund = [abundances[2]*(1.-fractions[3][0]-fractions[3][1]-fractions[3][2]),
             abundances[2]*fractions[3][0],
             abundances[3]*(1.-fractions[4][0]-fractions[4][1]),
             abundances[3]*fractions[4][0],
             abundances[3]*fractions[4][1],
             abundances[4]*fractions[5][1],
             abundances[5]*(1.-fractions[6][0]-fractions[6][1]-fractions[6][2]),
             abundances[5]*fractions[6][0],
             abundances[1]*(1.-fractions[2][0]-fractions[2][1]),
             abundances[1]*fractions[2][0],
             abundances[2]*fractions[3][1],
             abundances[5]*fractions[5][0]]

    linestrengths = lines.linestr(abund)

    t4 = temperature * 1.e-4
    emissivities["Hbeta"] = ne * nhp * 1.24e-38 * t4**-0.878
    emissivities["HeI_5876"] = ne * nhep * 1.69e-38 * t4**-1.065
    emissivities["OI_6300"] = density * linestrengths[0]

  return emissivities
