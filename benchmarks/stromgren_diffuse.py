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
# @file stromgren_diffuse.py
#
# @brief Read all stromgren_diffuse_*.hdf5 files in the directory and plot
# their hydrogen neutral fractions as a function of radius in a file
# stromgren_diffuse_*.png.
#
# We also plot the analytically calculated Stromgren radius for reference.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# load some libraries
import numpy as np
import h5py
import scipy.stats as stats
import glob
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True

# calculate the stromgren radius
# set the values for the parameters to the values used in the parameter file
alphaH = 4.0e-19  # m^3 s^-1
nH = 1.0e8  # m^-3
Q = 4.26e49  # s^-1

# this is the value for the reemission probability at 8,000 K
PR = 0.36392015

# adjust the total luminosity to take into account the extra luminosity due to
# the reemission
Q /= 1.0 - PR

# compute the Stromgren radius
Rs = (0.75 * Q / (np.pi * nH ** 2 * alphaH)) ** (1.0 / 3.0)

# output distances in pc
pc = 3.086e16  # m

Rs /= pc

# loop over all snapshot files in the directory
for fname in sorted(glob.glob("stromgren_diffuse_*.hdf5")):
    # give the user some clue about what is going on
    print("Processing", fname, "...")
    # open the snapshot file
    file = h5py.File(fname, "r")

    # read in the box size from the snapshot file
    box = np.array(file["/Header"].attrs["BoxSize"])
    box_center = 0.5 * box

    # read in the coordinates and neutral fractions
    coords = np.array(file["/PartType0/Coordinates"])
    nfracH = np.array(file["/PartType0/NeutralFractionH"])

    # calculate radii
    radius = np.sqrt(
        (coords[:, 0] - box_center[0]) ** 2
        + (coords[:, 1] - box_center[1]) ** 2
        + (coords[:, 2] - box_center[2]) ** 2
    )
    radius /= pc

    # bin the neutral fractions
    max_radius = max(radius)
    r_bin_edge = np.arange(0.0, max_radius, 0.02 * max_radius)
    r_bin = 0.5 * (r_bin_edge[1:] + r_bin_edge[:-1])
    nfracH_bin, _, _ = stats.binned_statistic(
        radius, nfracH, statistic="mean", bins=r_bin_edge
    )
    nfracH2_bin, _, _ = stats.binned_statistic(
        radius, nfracH ** 2, statistic="mean", bins=r_bin_edge
    )
    nfracH_sigma_bin = np.sqrt(nfracH2_bin - nfracH_bin ** 2)

    # plot the stromgren radius for reference
    pl.plot([Rs, Rs], [1.0e-7, 1.0], "r--", label=r"Str\"{o}mgren radius")
    # plot neutral fraction as a function of radius
    pl.semilogy(radius, nfracH, ".", color="grey", markersize=0.5, alpha=0.5)
    # plot the binned neutral fractions with error bar
    # we need to set the zorder to make sure the error bars are on top of the data
    # points
    pl.errorbar(
        r_bin,
        nfracH_bin,
        yerr=nfracH_sigma_bin,
        fmt=".",
        color="b",
        zorder=3,
        label=r"{\sc{}CMacIonize}",
    )

    # labels, legend and formatting...
    pl.xlabel(r"$r$ (pc)")
    pl.ylabel(r"$x_{\rm{}H}$")
    pl.legend(loc="best")
    pl.tight_layout()
    # save the file and reset the plot for the next file
    pl.savefig("{name}.png".format(name=fname[:-5]))
    pl.close()
