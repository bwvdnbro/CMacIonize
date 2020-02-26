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
# @file lexingtonHII40.py
#
# @brief This script reads in the last snapshot from the Lexington HII 40
# benchmark run and plots the neutral fractions, ionic fractions and temperature
# as a function of radius. The 6 panel plots are almost identical to those in
# Figure 1 of Wood, Mathis & Ercolano (2004), except that we replaced the
# temperature profile with sulphur abundance profiles and moved the temperature
# to a separate plot.
#
# We bin the results in 100 radial bins and plot the average with an error bar
# showing the standard deviation within each bin. The limits of the axes are
# set to the same limits as in the paper.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import some modules
import numpy as np
import h5py
import glob
import scipy.stats as stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

##
# @brief Bin the given quantity in the given radial bins.
#
# @param r Radii for the quantity to bin.
# @param q Quantity to bin.
# @param r_bin_edge Edges of the radial bins to use. Should contain one more
# value than the desired number of bins.
# @return Binned values and standard deviations.
##
def bin_quantity(r, q, r_bin_edge):
    q_bin, _, _ = stats.binned_statistic(
        r, q, statistic="mean", bins=r_bin_edge
    )
    q2_bin, _, _ = stats.binned_statistic(
        r, q ** 2, statistic="mean", bins=r_bin_edge
    )
    q_sigma_bin = np.sqrt(q2_bin - q_bin ** 2)
    return q_bin, q_sigma_bin


# define the parsec
pc = 3.086e16  # in m

# set the radial limits for the x axis
xlims = [1.0, 5.5]

# set the number of radial bins
numbin = 100

# get the last snapshot in the directory
f = sorted(glob.glob("lexingtonHII40_*.hdf5"))[-1]
print("Processing last snapshot:", f)

print("Reading file")
# open the file
file = h5py.File(f, "r")
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
ifracSp2 = np.array(file["/PartType0/NeutralFractionS+"])
ifracSp3 = np.array(file["/PartType0/NeutralFractionS++"])
ifracSp4 = np.array(file["/PartType0/NeutralFractionS+++"])
temp = np.array(file["/PartType0/Temperature"])

# compute additional data arrays
ifracCp1 = 1.0 - ifracCp2 - ifracCp3
ifracSp1 = 1.0 - ifracSp2 - ifracSp3 - ifracSp4

# read in the box size
box = np.array(file["/Header"].attrs["BoxSize"])

print("Done reading file")

print("Setting up radial bins")
# compute the radii
radius = np.sqrt(
    (coords[:, 0] - 0.5 * box[0]) ** 2
    + (coords[:, 1] - 0.5 * box[1]) ** 2
    + (coords[:, 2] - 0.5 * box[2]) ** 2
)
# convert them from m to pc
radius /= pc

# set up the radial bins...
rbin = np.linspace(xlims[0], xlims[1], numbin + 1)
rbin = np.concatenate(([0.0], rbin))
rbin = np.concatenate((rbin, [np.max(radius) * 1.1]))
# ...and the centres of the bins
rmid = 0.5 * (rbin[:-1] + rbin[1:])

print("Binning data")
# get the binned values and standard deviations
nfracHb, nfracHs = bin_quantity(radius, nfracH, rbin)
nfracHeb, nfracHes = bin_quantity(radius, nfracHe, rbin)
ifracOp1b, ifracOp1s = bin_quantity(radius, ifracOp1, rbin)
ifracOp2b, ifracOp2s = bin_quantity(radius, ifracOp2, rbin)
ifracCp1b, ifracCp1s = bin_quantity(radius, ifracCp1, rbin)
ifracCp2b, ifracCp2s = bin_quantity(radius, ifracCp2, rbin)
ifracCp3b, ifracCp3s = bin_quantity(radius, ifracCp3, rbin)
ifracNp1b, ifracNp1s = bin_quantity(radius, ifracNp1, rbin)
ifracNp2b, ifracNp2s = bin_quantity(radius, ifracNp2, rbin)
ifracNp3b, ifracNp3s = bin_quantity(radius, ifracNp3, rbin)
ifracNep1b, ifracNep1s = bin_quantity(radius, ifracNep1, rbin)
ifracNep2b, ifracNep2s = bin_quantity(radius, ifracNep2, rbin)
ifracSp1b, ifracSp1s = bin_quantity(radius, ifracSp1, rbin)
ifracSp2b, ifracSp2s = bin_quantity(radius, ifracSp2, rbin)
ifracSp3b, ifracSp3s = bin_quantity(radius, ifracSp3, rbin)
ifracSp4b, ifracSp4s = bin_quantity(radius, ifracSp4, rbin)
tempb, temps = bin_quantity(radius, temp, rbin)

print("Creating figures")
# create the abundances figure
fig, ax = pl.subplots(2, 3, figsize=(16, 12))

# we need to clip the errorbars that end up outside the figure box, since
# otherwise they are not shown properly
ax[0][0].set_yscale("log", nonposy="clip")
ax[0][0].errorbar(rmid, nfracHb, yerr=nfracHs, fmt=".", label="H0/H")
ax[0][0].errorbar(rmid, nfracHeb, yerr=nfracHes, fmt=".", label="He0/He")
ax[0][0].set_xlim(xlims[0], xlims[1])
ax[0][0].set_ylim(1.0e-5, 2.0)
ax[0][0].set_xlabel("$r$ (pc)")
ax[0][0].set_ylabel("Neutral fraction")
ax[0][0].legend(loc="best")

ax[0][1].set_yscale("log", nonposy="clip")
ax[0][1].errorbar(rmid, ifracOp1b, yerr=ifracOp1s, fmt=".", label="O+/O")
ax[0][1].errorbar(rmid, ifracOp2b, yerr=ifracOp2s, fmt=".", label="O++/O")
ax[0][1].set_xlim(xlims[0], xlims[1])
ax[0][1].set_ylim(1.0e-3, 2.0)
ax[0][1].set_xlabel("$r$ (pc)")
ax[0][1].set_ylabel("Ion fraction")
ax[0][1].legend(loc="best")

ax[0][2].set_yscale("log", nonposy="clip")
ax[0][2].errorbar(rmid, ifracCp1b, yerr=ifracCp1s, fmt=".", label="C+/C")
ax[0][2].errorbar(rmid, ifracCp2b, yerr=ifracCp2s, fmt=".", label="C++/C")
ax[0][2].errorbar(rmid, ifracCp3b, yerr=ifracCp3s, fmt=".", label="C+++/C")
ax[0][2].set_xlim(xlims[0], xlims[1])
ax[0][2].set_ylim(0.001, 2.0)
ax[0][2].set_xlabel("$r$ (pc)")
ax[0][2].set_ylabel("Ion fraction")
ax[0][2].legend(loc="best")

ax[1][0].set_yscale("log", nonposy="clip")
ax[1][0].errorbar(rmid, ifracNp1b, yerr=ifracNp1s, fmt=".", label="N+/N")
ax[1][0].errorbar(rmid, ifracNp2b, yerr=ifracNp2s, fmt=".", label="N++/N")
ax[1][0].errorbar(rmid, ifracNp3b, yerr=ifracNp3s, fmt=".", label="N+++/N")
ax[1][0].set_xlim(xlims[0], xlims[1])
ax[1][0].set_ylim(0.001, 2.0)
ax[1][0].set_xlabel("$r$ (pc)")
ax[1][0].set_ylabel("Ion fraction")
ax[1][0].legend(loc="best")

ax[1][1].set_yscale("log", nonposy="clip")
ax[1][1].errorbar(rmid, ifracNep1b, yerr=ifracNep1s, fmt=".", label="Ne+/Ne")
ax[1][1].errorbar(rmid, ifracNep2b, yerr=ifracNep2s, fmt=".", label="Ne++/Ne")
ax[1][1].set_xlim(xlims[0], xlims[1])
ax[1][1].set_ylim(0.001, 2.0)
ax[1][1].set_xlabel("$r$ (pc)")
ax[1][1].set_ylabel("Ion fraction")
ax[1][1].legend(loc="best")

ax[1][2].set_yscale("log", nonposy="clip")
ax[1][2].errorbar(rmid, ifracSp1b, yerr=ifracSp1s, fmt=".", label="S+/S")
ax[1][2].errorbar(rmid, ifracSp2b, yerr=ifracSp2s, fmt=".", label="S++/S")
ax[1][2].errorbar(rmid, ifracSp3b, yerr=ifracSp3s, fmt=".", label="S+++/S")
ax[1][2].errorbar(rmid, ifracSp4b, yerr=ifracSp4s, fmt=".", label="S++++/S")
ax[1][2].set_xlim(xlims[0], xlims[1])
ax[1][2].set_ylim(1.0e-6, 2.0)
ax[1][2].set_xlabel("$r$ (pc)")
ax[1][2].set_ylabel("Ion fraction")
ax[1][2].legend(loc="best")

# let matplotlib handle the whitespace in the figure to make it look nicer
pl.tight_layout()
# save the figure
abundances_filename = "{name}_abundances.png".format(name=f[:-5])
pl.savefig(abundances_filename)
print("Wrote", abundances_filename)
# close the plot, as we want to make a new one
pl.close()

# plot the temperature
pl.errorbar(rmid, tempb, yerr=temps, fmt=".")
pl.xlim(xlims[0], xlims[1])
pl.ylim(6000.0, 11000.0)
pl.ylabel("$T$ (K)")
pl.xlabel("$r$ (pc)")
pl.tight_layout()
temperature_filename = "{name}_temperature.png".format(name=f[:-5])
pl.savefig(temperature_filename)
print("Wrote", temperature_filename)
