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
# @file starbench.py
#
# @brief Read all starbench_*.hdf5 files in the directory and plot their density
# as a function of radius in a file starbench_*.png.
#
# We also plot the analytically calculated reference radii.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# load some libraries
import numpy as np
import h5py

# make sure matplotlib uses a non X backend for remote plotting
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import glob

# set some plot options
pl.rcParams["figure.figsize"] = (6, 8)
pl.rcParams["text.usetex"] = True

# set some conversion constants
pc = 3.086e16  # m
Myr = 3600.0 * 24.0 * 365.25 * 1.0e6  # s

# constant benchmark parameters
Rst = 0.314  # pc
ci = 12.85  # km/s

# convert to SI
Rst *= pc  # m
ci *= 1.0e3  # m/s

##
# @brief Spitzer reference solution.
#
# @param t Time (in s).
# @return Reference radius (in m).
##
def spitzer(t):
    global Rst, ci
    return Rst * (1.0 + 1.75 * ci * t / Rst) ** (4.0 / 7.0)


##
# @brief Hosokawa-Inutsuka reference solution.
#
# @param t Time (in s).
# @return Reference radius (in m).
##
def hosokawa_inutsuka(t):
    global Rst, ci
    return Rst * (1.0 + 1.75 * np.sqrt(4.0 / 3.0) * ci * t / Rst) ** (4.0 / 7.0)


# get the reference curves
trange = np.arange(0.0, 0.15 * Myr, 0.001 * Myr)
spitzerrange = np.array([spitzer(t) for t in trange])
hosinutrange = np.array([hosokawa_inutsuka(t) for t in trange])

# save the time and radius values for intermediate snapshot files
times = []
radii = []
# loop over all snapshot files in order
for f in sorted(glob.glob("starbench_*.hdf5")):
    # provide the user with some progress information
    print("plotting", f)
    # open the file and get the relevant data sets
    file = h5py.File(f, "r")
    coords = np.array(file["/PartType0/Coordinates"])
    rho = np.array(file["/PartType0/Density"])
    neutfracH = np.array(file["/PartType0/NeutralFractionH"])
    box = np.array(file["/Header"].attrs["BoxSize"])

    # compute the radii
    radius = np.sqrt(
        (coords[:, 0] - 0.5 * box[0]) ** 2
        + (coords[:, 1] - 0.5 * box[1]) ** 2
        + (coords[:, 2] - 0.5 * box[2]) ** 2
    )

    # get the time stamp from the snapshot
    times.append(file["/Header"].attrs["Time"] / Myr)

    # compute the average radius for cells with neutral fractions in the range
    # [0.8, 0.9]
    mean_transition_radius = 0.0
    num_transition = 0
    for i in range(len(neutfracH)):
        if neutfracH[i] < 0.9 and neutfracH[i] > 0.8:
            mean_transition_radius += radius[i]
            num_transition += 1
    if num_transition > 0:
        mean_transition_radius /= num_transition
    radii.append(mean_transition_radius / pc)

    # plot the density profile and ionization front position
    fig, ax = pl.subplots(2, 1)

    ax[0].plot(radius / pc, rho * 1.0e-3, "k.")
    ax[0].set_ylim(0.0, 3.0e-20)
    ax[0].set_xlabel("radius (pc)")
    ax[0].set_ylabel(r"density (g cm$^{-3}$)")

    ax[1].plot(trange / Myr, spitzerrange / pc, "r-")
    ax[1].plot(trange / Myr, hosinutrange / pc, "g-")
    ax[1].plot(times, radii, "kx")
    ax[1].set_ylim(0.0, 1.4)
    ax[1].plot(times[-1], radii[-1], "bx", markersize=10)
    ax[1].set_xlabel("time (Myr)")
    ax[1].set_ylabel("ionization front radius (pc)")

    pl.tight_layout()
    pl.savefig("{filename}.png".format(filename=f[:-5]))
    pl.close()

# plot the evolution of the ionization front
radii = np.array(radii)
times = np.array(times)

pl.plot(trange / Myr, spitzerrange / pc, "r-", label="Spitzer solution")
pl.plot(
    trange / Myr, hosinutrange / pc, "g-", label="Hosokawa-Inutsuka solution"
)
pl.plot(times[1:], radii[1:], "kx", label=r"{\sc{}CMacIonize}")
pl.xlabel("t (Myr)")
pl.ylabel("ionization front (pc)")
pl.legend(loc="best")
pl.tight_layout()
pl.savefig("starbench_profile.png")
