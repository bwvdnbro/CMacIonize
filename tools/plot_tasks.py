################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# import modules
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib.ticker import AutoMinorLocator
import pylab as pl
import sys
import argparse

# parse the command line arguments
argparser = argparse.ArgumentParser(
    description="Plot task plot based on a given task output file."
)

argparser.add_argument("-n", "--name", action="store", required=True)
argparser.add_argument("-l", "--labels", action="store_true")

args = argparser.parse_args(sys.argv[1:])

name = args.name

# change the default matplotlib settings to get nicer plots
pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (12, 10)
pl.rcParams["font.size"] = 14

# name labels and colours for the various task types
# add extra colours and labels here if new tasks are created
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
task_colors = pl.cm.ScalarMappable(cmap="tab20").to_rgba(
    np.linspace(0.0, 1.0, len(task_names))
)

# load the program time data
try:
    ptime = np.loadtxt("program_time.txt")
    if len(ptime.shape) > 1:
        ptime = np.array(sorted(ptime, key=lambda line: line[0]))
    else:
        ptime = np.array([ptime])
except:
    ptime = np.array([[0, -1, -1, 1.0]])

# load the data
print("Plotting tasks for", name, "...")
data = np.loadtxt(name)

if ptime[0, 1] < 0:
    ptime[0, 1] = data[:, 2].min()
    ptime[0, 2] = data[:, 3].max()
if (data[:, 4] == -1).sum() == 0:
    data = np.append(data, [[0, 0, ptime[0, 1], ptime[0, 2], -1, 0]], axis=0)

task_flags = [
    len(data[data[:, 4] == task]) > 0 for task in range(len(task_names))
]

# get information about the system
nthread = int(data[:, 1].max()) + 1
nproc = int(data[:, 0].max()) + 1

# get the minimum and maximum time stamp and compute the time to fraction
# conversion factor for each node
tmin = np.zeros(nproc)
tmin_in_s = np.zeros(nproc)
tmax = np.zeros(nproc)
tmax_in_s = np.zeros(nproc)
tconv = np.zeros(nproc)
tconv_in_s = np.zeros(nproc)
for iproc in range(nproc):
    procline = data[(data[:, 0] == iproc) & (data[:, 4] == -1)]
    if len(procline) > 1:
        print("Too many node information lines!")
        exit()
    tmin[iproc] = procline[0, 2]
    tmin_in_s[iproc] = (
        (tmin[iproc] - ptime[iproc][1])
        * ptime[iproc][3]
        / (ptime[iproc][2] - ptime[iproc][1])
    )
    tmax[iproc] = procline[0, 3]
    tmax_in_s[iproc] = (
        (tmax[iproc] - ptime[iproc][1])
        * ptime[iproc][3]
        / (ptime[iproc][2] - ptime[iproc][1])
    )
    tconv[iproc] = 1.0 / (tmax[iproc] - tmin[iproc])
    tconv_in_s[iproc] = (
        (tmax[iproc] - tmin[iproc])
        * ptime[iproc][3]
        / (ptime[iproc][2] - ptime[iproc][1])
    )

ttot_min = tmin_in_s.min()
ttot_max = tmax_in_s.max()

## make the plot

fig, ax = pl.subplots(1, 1, sharex=True)

ax.set_xlim(
    ttot_min - 0.05 * (ttot_max - ttot_min),
    ttot_max + 0.05 * (ttot_max - ttot_min),
)
ax.axvline(x=ttot_min, linestyle="--", color="k", linewidth=0.8)
ax.axvline(x=ttot_max, linestyle="--", color="k", linewidth=0.8)

# now plot the tasks
alltime = 0
# loop over the processes
for iproc in range(nproc):
    # filter out the data for this process
    process = data[(data[:, 0] == iproc) & (data[:, 4] != -1)]

    # loop over the threads
    for i in range(nthread):
        # filter out the data for this thread
        thread = process[(process[:, 1] == i)][:, 1:]

        # create the task plot
        bar = [
            (
                (task[1] - tmin[iproc]) * tconv[iproc],
                (task[2] - task[1]) * tconv[iproc],
            )
            for task in thread
        ]
        if len(thread[0]) > 4:
            tottime = thread[:, 4].sum() * tconv[iproc]
        else:
            tottime = np.array([line[1] for line in bar]).sum()
        alltime += tottime
        bar = [
            (
                line[0] * tconv_in_s[iproc] + tmin_in_s[iproc],
                line[1] * tconv_in_s[iproc],
            )
            for line in bar
        ]
        colors = [task_colors[int(task[3])] for task in thread]
        ax.broken_barh(
            bar,
            (iproc * nthread + i - 0.4, 0.8),
            facecolors=colors,
            edgecolor="none",
        )
        # optionally add labels
        if args.labels:
            # status text
            label = ""
            if nproc > 1:
                label += "rank {0} - ".format(iproc)
            if nthread > 1:
                label += "thread {0} - ".format(i)
            label += "{0:.2f} \% load".format(tottime * 100.0)
            ax.text(
                0.5 * (ttot_min + ttot_max),
                iproc * nthread + i + 0.2,
                label,
                ha="center",
                bbox=dict(facecolor="white", alpha=0.9),
            )
            # per task fraction text
            label = ""
            for itask in range(len(task_colors)):
                if task_flags[itask]:
                    tottime = np.array(
                        [
                            (task[2] - task[1]) * tconv[iproc]
                            for task in thread
                            if task[3] == itask
                        ]
                    ).sum()
                    label += "{0}: {1:.2f} \% - ".format(
                        task_names[itask], tottime * 100.0
                    )
            ax.text(
                0.5 * (ttot_min + ttot_max),
                iproc * nthread + i - 0.2,
                label[:-2],
                ha="center",
                bbox=dict(facecolor="white", alpha=0.9),
            )

# add empty blocks for the legend
for i in range(len(task_colors)):
    if task_flags[i]:
        ax.plot([], [], color=task_colors[i], label=task_names[i])

# add the legend and clean up the axes
ax.legend(loc="upper center", ncol=min(3, len(task_colors)))
ax.set_ylim(-1.0, nproc * nthread * 1.1)
ax.set_yticks([])

# get the total idle time and add information to the title
alltime /= nthread * nproc
ax.set_title("Total empty fraction: {0:.2f} \%".format((1.0 - alltime) * 100.0))
ax.set_xlabel("Simulation time (s)")

ax.xaxis.set_minor_locator(AutoMinorLocator())

for iproc in range(1, nproc):
    ax.axhline(y=iproc * nthread - 0.5, linestyle="-", linewidth=0.8, color="k")

# finalize and save the plot
pl.tight_layout()
pl.savefig("{0}.png".format(name[:-4]))
