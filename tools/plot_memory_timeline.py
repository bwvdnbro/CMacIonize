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
# @file plot_memory_timeline.py
#
# @brief Script to plot a time line of the memory usage of a task based parallel
# run.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules
import numpy as np
import matplotlib

matplotlib.use("Agg")
import pylab as pl

# use TeX in plots
pl.rcParams["text.usetex"] = True

# pre-define the MB
MB = 1 << 20

# load the memory usage time line data
data = np.loadtxt(
    "memory_timeline.txt",
    delimiter="\t",
    dtype={
        "names": ("label", "virtual size", "physical size", "timestamp"),
        "formats": ("S100", "u8", "u8", "u8"),
    },
)

# convert sizes to MB
for row in data:
    row["virtual size"] /= MB
    row["physical size"] /= MB

# plot the time line
t = np.arange(0, len(data))
pl.plot(t, data["virtual size"], label="virtual")
pl.plot(t, data["physical size"], label="physical")

# appropriate labels
pl.ylabel("memory usage (MB)")

# use snapshot labels as x-label
pl.gca().set_xticks(t)
pl.gca().set_xticklabels(data["label"], rotation="vertical")

pl.legend(loc="best")

# finalise and save image
pl.tight_layout()
pl.savefig("memory_timeline.png", dpi=300, bbox_inches="tight")
