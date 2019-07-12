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
# @file plot_tasks.py
#
# @brief Script to plot the task plot for a given file with task output.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import matplotlib

matplotlib.use("Agg")
import pylab as pl
import argparse

# parse the command line arguments
argparser = argparse.ArgumentParser(
    description="Plot task time stats based on a given task output file."
)

argparser.add_argument("-n", "--name", action="store", required=True)
argparser.add_argument("-o", "--output", action="store", required=True)
argparser.add_argument("-m", "--max", action="store", default=-1, type=float)
argparser.add_argument("-b", "--bins", action="store", default=100, type=int)
argparser.add_argument("-y", "--ylim", action="store", default=-1, type=float)

args = argparser.parse_args()

# list of all supported tasks
task_names = [
    "source photon (discrete)",
    "source photon (continuous)",
    "photon traversal",
    "reemission",
    "temperature/ionization state",
    "send",
    "receive",
    "gradsweep internal",
    "gradsweep neighbour",
    "gradsweep boundary",
    "slope limiter",
    "predict primitives",
    "fluxsweep internal",
    "fluxsweep neighbour",
    "fluxsweep boundary",
    "update conserved",
    "update primitives",
]

# load the task data
data = np.loadtxt(args.name)
tasks = []
labels = []
taskmax = 0
# loop over all task types
for i in range(len(task_names)):
    # filter out all tasks of this type
    task = data[data[:, 4] == i]

    # skip tasks that are not present
    if len(task) == 0:
        continue

    # compute the time spent in each task
    times = task[:, 3] - task[:, 2]
    taskmax = max(taskmax, times.max())

    # add to the lists
    tasks.append(times)
    labels.append(task_names[i])

# plot a histogram of all tasks that are present
if args.max > 0:
    taskmax = args.max
pl.hist(tasks, np.linspace(0, taskmax, args.bins + 1), label=labels)

if args.ylim > 0:
    pl.ylim(0, args.ylim)

pl.legend(loc="best")
pl.savefig(args.output, dpi=300)
