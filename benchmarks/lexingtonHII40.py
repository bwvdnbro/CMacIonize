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
# @file lexingtonHII40.py
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import h5py
import pylab as pl
import glob

for f in sorted(glob.glob("lexingtonHII40*.hdf5")):
  print "processing", f, "..."

  file = h5py.File(f, 'r')
  coords = np.array(file["/PartType0/Coordinates"])
  nfracH = np.array(file["/PartType0/NeutralFractionH"])

  box = np.array(file["/Header"].attrs["BoxSize"])

  radius = np.sqrt((coords[:,0] - 0.5*box[0])**2 +
                   (coords[:,1] - 0.5*box[1])**2 +
                   (coords[:,2] - 0.5*box[2])**2)

  pl.semilogy(radius, nfracH, "kx")
  pl.savefig("{name}.png".format(name = f[:-5]))
  pl.close()
