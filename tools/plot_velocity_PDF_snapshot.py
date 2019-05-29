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
# @file plot_velocity_PDF_snapshot.py
#
# @brief Script that reads a velocity PDF snapshot file.
#
# The script takes two command line arguments: the name of the file to plot,
# and the name of the desired output file (should be an image format supported
# by Matplotlib). All other relevant information is automatically read from the
# output file.
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
args = argparser.parse_args()

# open the file and read the PDF information on the first 2 lines
file = open(args.file, "r")
min_max = np.fromstring(file.readline(), sep="\t", count=2)
bin_info = np.fromstring(file.readline(), sep="\t", count=1)
file.close()

# now parse the rest of the file into a single array
hist = np.loadtxt(args.file, skiprows=2)

number_of_bins = len(hist)
bin_edges = np.linspace(
    0.0, (number_of_bins + 1) * bin_info[0], number_of_bins + 1
)
bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])

pl.plot(bin_centres, hist)
pl.xlabel("$v$ (m s$^{-1}$)")
pl.ylabel("bin count")
pl.title("Velocity PDF (min: {0:.2e}, max: {1:.2e})".format(*min_max))

# save the image
pl.tight_layout()
pl.savefig(args.output, dpi=300, bbox_inches="tight")
