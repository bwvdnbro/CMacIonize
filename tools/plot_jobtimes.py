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
# @file get_parameterfile.py
#
# @brief Script that plots a bar plot of tasks executed on different threads in
# a shared memory parallel context, which is useful to identify serial parts of
# the code and load imbalances between different threads.
#
# The script takes a single command line parameter: the number of threads that
# should be plotted. The script then assumes that an equal number of
# corresponding 'jobtimes_X.txt' files exists in the current working directory,
# with X in [0, number of threads[.
#
# To produce these files, configure the code with the
# -DACTIVATE_OUTPUT_CYCLES=True option. Note that this has a considerable impact
# on the runtime and output size of the code, so this is only suitable for
# relatively small simulations.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules: pylab is used to make plots, sys is needed to parse the
# command line arguments
import pylab as pl
import sys

# global list of colours for the different job names reported by the code
jobs = {
    "testjob": "r",
    "densitygrid_traversal<N11DensityGrid33DensityGridInitializationFunctionE>": "b",
    "photonshootjob": "y",
    "densitygrid_traversal"
    "<N25IonizationStateCalculator33IonizationStateCalculatorFunctionE>": "g",
    "densitygrid_traversal"
    "<N25TemperatureCalculator33TemperatureCalculatorFunctionE>": "k",
    "fractaldensitymask_construction": "c",
    "voronoigrid_construction": "m",
}

##
# @brief Auxiliary function that parses a 'jobtimes_X.txt' file and returns a
# list of jobs and corresponding start and end cycle.
#
# @param filename Name to read.
# @return List of jobs and corresponding start and end cycle.
##
def get_times(filename):
    # open the file for reading (this displays a pretty clear error if the file
    # does not exist, so no need to check this ourselves)
    file = open(filename, "r")
    # read the lines from the file
    lines = file.readlines()
    # parse the lines and add them to a list
    times = []
    for line in lines:
        data = line.split()
        times.append([data[0], int(data[1]), int(data[2])])
    return times


##
# @brief Main script routine.
#
# @param args Command line options passed on to the script.
##
def main(args):
    # tell the interpreter we want to use the global 'jobs' list
    global jobs
    # tell the interpreter we want to use the 'pl' module
    global pl

    # parse the command line arguments
    # note that the first command line argument is always the name of the script
    if len(args) < 2:
        print("Usage: python plot_jobtimes.py NUMBER_OF_THREADS")
        exit()
    num_threads = int(args[1])

    # do the actual plotting
    # loop over the number of threads
    for i in range(num_threads):
        # get the data
        times = get_times("jobtimes_{i}.txt".format(i=i))
        # plot each data group in the corresponding colour
        for time in times:
            for key in jobs:
                if time[0] == key:
                    pl.barh(
                        i,
                        time[2] - time[1],
                        left=time[1] - times[0][1],
                        color=jobs[key],
                    )
    # show the plot
    # the program will resume when the window is closed by the user
    pl.show()


# make sure the main() routine is called
if __name__ == "__main__":
    main(sys.argv)
