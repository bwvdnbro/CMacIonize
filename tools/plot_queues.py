################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file plot_queues.py
#
# @brief Script to plot an overview of the queue load per iteration for
# the different processes and threads.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules
import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import glob
from operator import itemgetter

# get a list of all files that are present
files = sorted(glob.glob("queues_??.txt"))
# collect the data
alldata = []
for file in files:
  data = np.loadtxt(file)
  data = data[data[:,1] != -1]
  if len(data.shape) > 1:
    data = np.array(sorted(data, key = itemgetter(0, 1)))
  else:
    data = [data]
  alldata.append(data)

alldata = np.array(alldata)

# prepare data for plotting
nproc = int(alldata[:,:,0].max()) + 2
nthread = int(alldata[:,:,1].max()) + 1
alldata[:,:,0] = alldata[:,:,0] / nproc
alldata[:,:,1] = alldata[:,:,1] / (nproc * nthread)
steps = np.arange(len(files))

# plot the contributions for the different processes
for i in range(nproc - 1):
  for j in range(nthread):
    pl.bar(alldata[:,i * nthread + j,0] + alldata[:,i * nthread + j,1] + steps,
           alldata[:,i * nthread + j,2], 1. / (nproc * nthread))

# plot layout
pl.xticks([])
pl.ylabel("Number of used queue elements")

# save the plot
pl.tight_layout()
pl.savefig("queue_stats.png")
