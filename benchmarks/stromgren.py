#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @brief Read all stromgren_*.hdf5 files in the directory and plot their
# hydrogen neutral fractions as a function of radius in a file stromgren_*.png.
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
sigmaH = 6.3e-22  # m^2

# compute the Stromgren radius
Rs = (0.75 * Q / (np.pi * nH ** 2 * alphaH)) ** (1.0 / 3.0)

# compute the reference neutral fraction profile
rref = np.linspace(0.0, 1.2 * Rs, 1200)
xref = np.zeros(rref.shape)
integral = 0.0
factor = 0.125 * Q * sigmaH / (np.pi * nH * alphaH)
intfac = 0.0005 * Rs * nH * sigmaH
for i in range(1, len(rref)):
    A = factor * np.exp(-integral) / rref[i] ** 2
    xref[i] = 1.0 + A - np.sqrt(2.0 * A + A ** 2)
    integral += intfac * (xref[i - 1] + xref[i])

# output distances in pc
pc = 3.086e16  # m
Rs /= pc

# loop over all snapshot files in the directory
for fname in sorted(glob.glob("stromgren_*.hdf5")):
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
    pl.gca().axvline(
        x=Rs, color="r", linestyle="--", label=r"Str\"{o}mgren radius"
    )
    # plot neutral fraction as a function of radius
    pl.semilogy(radius, nfracH, ".", color="grey", markersize=0.5, alpha=0.5)
    # plot the reference curve
    pl.plot(rref / pc, xref, "r-", label="Reference profile")
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
