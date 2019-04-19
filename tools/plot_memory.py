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
# @file plot_memory.py
#
# @brief Script to plot a pie chart of the memory usage for a task based
# parallel run.
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

# predefine the MB (in bytes)
MB = 1<<20

# load the memory allocation data
data = np.loadtxt("memory.txt", delimiter = "\t",
                  dtype = {"names": ("label", "virtual size", "physical size",
                                     "timestamp"),
                           "formats": ("S100", "u8", "u8", "u8")})

# filter out small contributions (less than 1% of total memory usage)
memsum = data["virtual size"].sum()
data = data[data["virtual size"] > 0.01 * memsum]

# convert sizes to MB
for row in data:
  row["label"] += " ({0:.0f} MB)".format(float(row["virtual size"]) / MB)

# create the pie chart
pl.pie(data["virtual size"], explode = np.ones(len(data)) * 0.1,
       labels = data["label"])

# put the total memory usage in the title
pl.title("Total memory: {0:.0f} MB".format(memsum / MB))

# force aspect ratio
pl.gca().set_aspect("equal")

# save the image
pl.savefig("memory.png", dpi = 300, bbox_inches = "tight")
