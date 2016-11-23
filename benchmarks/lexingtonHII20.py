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
# @file lexingtonHII20.py
#
# @brief This script reads in the snapshots from the Lexington HII 20
# benchmark run and plots the neutral fractions, ionic fractions and temperature
# as a function of radius. The 6 panel plots are identical to those in Figure 2
# of Wood, Mathis & Ercolano (2004).
#
# We bin the results in 100 radial bins and plot the average with an error bar
# showing the standard deviation within each bin. The limits of the axes are
# set to the same limits as in the paper.
#
# We do not plot C+/C, since C+ is not included in our version of the code.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import some modules
import numpy as np
import h5py
import pylab as pl
import glob

##
# @brief Get the average values in the radial bins, and the standard deviations
# within each bin.
#
# @param values Array to bin.
# @param ibins Indices of the bin containing each value: values[i] is in bin
# ibins[i]. These should be precomputed using e.g. numpy.digitize.
# @param size Number of bins in total. This should be equal to the maximum index
# value present in ibins.
# @return Binned values and standard deviations.
##
def get_averages(values, ibins, size):
  bins = np.zeros(size)
  nbins = np.zeros(size)
  sbins = np.zeros(size)

  for i in range(len(values)):
    bins[ibins[i]-1] += values[i]
    nbins[ibins[i]-1] += 1
  bins /= nbins

  for i in range(len(values)):
    dv = values[i] - bins[ibins[i]-1]
    sbins[ibins[i]-1] += dv**2
  sbins /= nbins
  sbins = np.sqrt(sbins)

  return bins, sbins

# define the parsec
pc = 3.086e16 # in m

# set the radial limits for the x axis
xlims = [1., 3.5]

# set the number of radial bins
numbin = 100

# loop over all snapshots in the folder
for f in sorted(glob.glob("lexingtonHII20_*.hdf5")):
  # tell the user what we are doing
  print "processing", f, "..."

  # open the file
  file = h5py.File(f, 'r')
  # read in the data arrays
  coords = np.array(file["/PartType0/Coordinates"])
  nfracH = np.array(file["/PartType0/NeutralFractionH"])
  nfracHe = np.array(file["/PartType0/NeutralFractionHe"])
  ifracOp1 = np.array(file["/PartType0/NeutralFractionO"])
  ifracOp2 = np.array(file["/PartType0/NeutralFractionO+"])
  ifracCp2 = np.array(file["/PartType0/NeutralFractionC+"])
  ifracCp3 = np.array(file["/PartType0/NeutralFractionC++"])
  ifracNp1 = np.array(file["/PartType0/NeutralFractionN"])
  ifracNp2 = np.array(file["/PartType0/NeutralFractionN+"])
  ifracNp3 = np.array(file["/PartType0/NeutralFractionN++"])
  ifracNep1 = np.array(file["/PartType0/NeutralFractionNe"])
  ifracNep2 = np.array(file["/PartType0/NeutralFractionNe+"])
  temp = np.array(file["/PartType0/Temperature"])

  # read in the box size
  box = np.array(file["/Header"].attrs["BoxSize"])

  # compute the radii
  radius = np.sqrt((coords[:,0] - 0.5*box[0])**2 +
                   (coords[:,1] - 0.5*box[1])**2 +
                   (coords[:,2] - 0.5*box[2])**2)
  # convert them from m to pc
  radius /= pc

  # set up the radial bins...
  rbin = np.linspace(xlims[0], xlims[1], numbin+1)
  rbin = np.concatenate(([0.], rbin))
  rbin = np.concatenate((rbin, [np.max(radius)*1.1]))
  # ...and the centres of the bins
  rmid = 0.5*(rbin[:-1] + rbin[1:])

  # precalculate the bin indices
  ibins = np.digitize(radius, rbin, right = False)

  # get the binned values and standard deviations
  nfracHb, nfracHs = get_averages(nfracH, ibins, numbin+2)
  nfracHeb, nfracHes = get_averages(nfracHe, ibins, numbin+2)
  ifracOp1b, ifracOp1s = get_averages(ifracOp1, ibins, numbin+2)
  ifracOp2b, ifracOp2s = get_averages(ifracOp2, ibins, numbin+2)
  ifracCp2b, ifracCp2s = get_averages(ifracCp2, ibins, numbin+2)
  ifracCp3b, ifracCp3s = get_averages(ifracCp3, ibins, numbin+2)
  ifracNp1b, ifracNp1s = get_averages(ifracNp1, ibins, numbin+2)
  ifracNp2b, ifracNp2s = get_averages(ifracNp2, ibins, numbin+2)
  ifracNp3b, ifracNp3s = get_averages(ifracNp3, ibins, numbin+2)
  ifracNep1b, ifracNep1s = get_averages(ifracNep1, ibins, numbin+2)
  ifracNep2b, ifracNep2s = get_averages(ifracNep2, ibins, numbin+2)
  tempb, temps = get_averages(temp, ibins, numbin+2)

  # create the figure
  fig, ax = pl.subplots(2, 3, figsize = (16, 12))

  ax[0][0].set_yscale("log", nonposy = "clip")
  ax[0][0].errorbar(rmid, nfracHb, yerr = nfracHs, fmt = "o", label = "H0/H")
  ax[0][0].errorbar(rmid, nfracHeb, yerr = nfracHes, fmt = "o",
                    label = "He0/He")
  ax[0][0].set_xlim(xlims[0], xlims[1])
  ax[0][0].set_ylim(1.e-4, 2.)
  ax[0][0].set_xlabel("$r$ (pc)")
  ax[0][0].set_ylabel("Neutral fraction")
  ax[0][0].legend(loc = "best")

  ax[0][1].set_yscale("log", nonposy = "clip")
  ax[0][1].errorbar(rmid, ifracOp1b, yerr = ifracOp1s, fmt = "o",
                    label = "O+/O")
  ax[0][1].errorbar(rmid, ifracOp2b, yerr = ifracOp2s, fmt = "o",
                    label = "O++/O")
  ax[0][1].set_xlim(xlims[0], xlims[1])
  ax[0][1].set_ylim(1.e-3, 2.)
  ax[0][1].set_xlabel("$r$ (pc)")
  ax[0][1].set_ylabel("Ion fraction")
  ax[0][1].legend(loc = "best")

  ax[0][2].set_yscale("log", nonposy = "clip")
  ax[0][2].errorbar(rmid, ifracCp2b, yerr = ifracCp2s, fmt = "o",
                    label = "C++/C")
  ax[0][2].errorbar(rmid, ifracCp3b, yerr = ifracCp3s, fmt = "o",
                    label = "C+++/C")
  ax[0][2].set_xlim(xlims[0], xlims[1])
  ax[0][2].set_ylim(0.001, 2.)
  ax[0][2].set_xlabel("$r$ (pc)")
  ax[0][2].set_ylabel("Ion fraction")
  ax[0][2].legend(loc = "best")

  ax[1][0].set_yscale("log", nonposy = "clip")
  ax[1][0].errorbar(rmid, ifracNp1b, yerr = ifracNp1s, fmt = "o",
                    label = "N+/N")
  ax[1][0].errorbar(rmid, ifracNp2b, yerr = ifracNp2s, fmt = "o",
                    label = "N++/N")
  ax[1][0].errorbar(rmid, ifracNp3b, yerr = ifracNp3s, fmt = "o",
                    label = "N+++/N")
  ax[1][0].set_xlim(xlims[0], xlims[1])
  ax[1][0].set_ylim(0.001, 2.)
  ax[1][0].set_xlabel("$r$ (pc)")
  ax[1][0].set_ylabel("Ion fraction")
  ax[1][0].legend(loc = "best")

  ax[1][1].set_yscale("log", nonposy = "clip")
  ax[1][1].errorbar(rmid, ifracNep1b, yerr = ifracNep1s, fmt = "o",
                    label = "Ne+/Ne")
  ax[1][1].errorbar(rmid, ifracNep2b, yerr = ifracNep2s, fmt = "o",
                    label = "Ne++/Ne")
  ax[1][1].set_xlim(xlims[0], xlims[1])
  ax[1][1].set_ylim(0.001, 2.)
  ax[1][1].set_xlabel("$r$ (pc)")
  ax[1][1].set_ylabel("Ion fraction")
  ax[1][1].legend(loc = "best")

  ax[1][2].errorbar(rmid, tempb, yerr = temps, fmt = "o")
  ax[1][2].set_xlim(xlims[0], xlims[1])
  ax[1][2].set_ylim(5000., 9500.)
  ax[1][2].set_xlabel("$r$ (pc)")
  ax[1][2].set_ylabel("Temperature (K)")

  pl.tight_layout()
  pl.savefig("{name}.png".format(name = f[:-5]))
  pl.close()
