################################################################################
# This file is part of CMacIonize
# Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file plot_surface_density_snapshot.py
#
# @brief Script that generates an image for a surface density output file.
#
# The script takes two command line arguments: the name of the file to plot,
# and the name of the desired output file (should be an image format supported
# by Matplotlib). All other relevant information is automatically read from the
# output file.
# An additional optional command line argument (--log) can be used to plot
# the surface density on a logarithmic rather than a linear scale.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules:
#  - numpy for file reading and array operations
#  - matplotlib for plotting (using the Agg backend)
#  - argparse for smart command line argument parsing
import numpy as np
import matplotlib

matplotlib.use("Agg")
import pylab as pl
import argparse

# enable LaTeX in the matplotlib output
pl.rcParams["text.usetex"] = True

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("--file", "-f", action="store", required=True)
argparser.add_argument("--output", "-o", action="store", required=True)
argparser.add_argument("--log", "-l", action="store_true")
args = argparser.parse_args()

# open the file and read the image layout information on the first 4 lines
file = open(args.file, "r")
ngrid = np.fromstring(file.readline(), sep="\t", count=2, dtype=int)
ncell = np.fromstring(file.readline(), sep="\t", count=2, dtype=int)
box_anchor = np.fromstring(file.readline(), sep="\t", count=2)
box_sides = np.fromstring(file.readline(), sep="\t", count=2)
file.close()

# now parse the rest of the file into a single array
pixels = np.loadtxt(args.file, skiprows=4)

# reorder the 1D array into a 2D image
blocksize = ncell[0] * ncell[1]
ntot = np.array([ngrid[0] * ncell[0], ngrid[1] * ncell[1]])
image = np.zeros((ntot[0], ntot[1]))
for ix in range(ngrid[0]):
    for iy in range(ngrid[1]):
        iblock = ix * ngrid[1] + iy
        image[
            ix * ncell[0] : (ix + 1) * ncell[0],
            iy * ncell[1] : (iy + 1) * ncell[1],
        ] = pixels[iblock * blocksize : (iblock + 1) * blocksize].reshape(
            (ncell[0], ncell[1])
        )

# set up the coordinate axes based on the box information in the snapshot
x = np.linspace(box_anchor[0], box_sides[0], ntot[0])
y = np.linspace(box_anchor[1], box_sides[1], ntot[1])

# plot the image
fig, ax = pl.subplots(1, 1)

if args.log:
    image = np.log10(image)
rhoplot = ax.contourf(x, y, image, 500)
ax.set_xlabel("$x$ (m)")
ax.set_ylabel("$y$ (m)")

# add a colorbar and overwrite default matplotlib ticks (and tick labels for
# logarithmic plots)
cbar = fig.colorbar(rhoplot, ax=ax, label="$\\rho{}$ (kg m$^{-3}$)")
rhomin = np.ceil(image.min())
rhomax = np.floor(image.max())
ticks = np.arange(rhomin, rhomax + 0.5, 1.0)
cbar.set_ticks(ticks)
if args.log:
    ticklabels = []
    for tick in ticks:
        ticklabels.append("$10^{{{0:.0f}}}$".format(tick))
    cbar.set_ticklabels(ticklabels)

# save the image
pl.tight_layout()
pl.savefig(args.output, dpi=300, bbox_inches="tight")
