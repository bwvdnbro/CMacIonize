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
# @file plot_timeline.py
#
# @brief Script that generates an image of the simulation time line.
#
# The script takes two required command line arguments: the name of the file to
# plot, and the name of the desired output file (should be an image format
# supported by Matplotlib). All other relevant information is automatically
# read from the output file.
# Additional optional command line arguments can be used to control the depth
# in the time log that is plotted (--min-depth, --max-depth).
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
argparser.add_argument("--min-depth", "-m", action="store", default=0, type=int)
argparser.add_argument(
    "--max-depth", "-M", action="store", default=-1, type=int
)
args = argparser.parse_args()

# load the time log
data = np.loadtxt(
    args.file,
    delimiter="\t",
    dtype={
        "names": ("id", "pid", "depth", "tic", "toc", "start", "end", "label"),
        "formats": ("u4", "u4", "u4", "u8", "u8", "f8", "f8", "S100"),
    },
)

# get the depth values
mindepth = args.min_depth
if args.max_depth < 0:
    maxdepth = data["depth"].max()

# plot all requested depths
for depth in range(mindepth, maxdepth + 1):
    # select data for this depth
    idx = data["depth"] == depth

    # create a bar plot
    bar = [(line["start"], line["end"] - line["start"]) for line in data[idx]]
    colors = ["C{0}".format(i % 10) for i in range(len(data[idx]))]
    pl.broken_barh(bar, (depth - 0.4, 0.8), facecolors=colors, edgecolor="none")
    # add labels
    labels = [line["label"] for line in data[idx]]
    for i in range(len(labels)):
        pl.text(
            bar[i][0] + 0.5 * bar[i][1],
            depth - 0.2,
            labels[i],
            ha="center",
            bbox=dict(facecolor="white", alpha=0.9),
        )

# disable vertical ticks
pl.gca().set_yticks([])
# set the x axis label
pl.xlabel("run time (s)")

# save the plot
pl.tight_layout()
pl.savefig(args.output, dpi=300, bbox_inches="tight")
