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
# @file stromgren.py
#
# @brief Read all snapshot*.hdf5 files in the directory and plot their hydrogen
# neutral fractions as a function of radius in a file snapshot*.png.
#
# We also plot the analytically calculated stromgren radius for reference.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# load some libraries
import numpy as np
import h5py
import pylab as pl
import glob

# calculate the stromgren radius
# we assume a temperature of 8000 K, a density of 100 cm^-3 and an ionizing
# luminosity of 4.26e49 s^-1.
alphaH = 4.91452e-19 # m^3 s^-1
nH = 1.e8 # m^-3
Q = 4.26e49 # s^-1
Rs = (0.75*Q/np.pi/nH**2/alphaH)**(1./3.)

# output distances in pc
pc = 3.086e16 # m

Rs /= pc

# loop over all snapshot files in the directory
for fname in sorted(glob.glob("snapshot*.hdf5")):
  # give the user some clue about what is going on
  print "Processing", fname, "..."
  # open the snapshot file
  file = h5py.File(fname, 'r')

  # read in the box size from the snapshot file
  box = np.array(file["/Header"].attrs["BoxSize"])
  box_center = 0.5*box

  # read in the coordinates and neutral fractions
  coords = np.array(file["/PartType0/Coordinates"])
  nfracH = np.array(file["/PartType0/NeutralFractionH"])

  # calculate radii
  radius = np.sqrt((coords[:,0]-box_center[0])**2 + 
                   (coords[:,1]-box_center[1])**2 +
                   (coords[:,2]-box_center[2])**2)
  radius /= pc

  # plot the stromgren radius for reference
  pl.plot([Rs, Rs], [1.e-7, 1.], "r--", label = "Stromgren radius")
  # plot neutral fraction as a function of radius
  pl.semilogy(radius, nfracH, "k.")
  
  # labels, legend and formatting...
  pl.xlabel(r"$r$ (pc)")
  pl.ylabel(r"$x_{\rm{}H}$")
  pl.legend(loc = "best")
  pl.tight_layout()
  # save the file and reset the plot for the next file
  pl.savefig("{name}.png".format(name = fname[:-5]))
  pl.close()
