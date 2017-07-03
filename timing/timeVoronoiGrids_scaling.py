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
# @file timeVoronoiGrids_scaling.py
#
# @brief Script to plot the scaling results of timeVoronoiGrids.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import pylab as pl

rand_old = np.loadtxt("timeVoronoiGrids_scaling_random_old.txt")
rand_new = np.loadtxt("timeVoronoiGrids_scaling_random_new.txt")
reg_old = np.loadtxt("timeVoronoiGrids_scaling_regular_old.txt")
reg_new = np.loadtxt("timeVoronoiGrids_scaling_regular_new.txt")

rand_old_speedup = rand_old[0,1]/rand_old[:,1]
rand_new_speedup = rand_new[0,1]/rand_new[:,1]
reg_old_speedup = reg_old[0,1]/reg_old[:,1]
reg_new_speedup = reg_new[0,1]/reg_new[:,1]
rand_old_speedup_std = rand_old_speedup * \
  np.sqrt((rand_old[0,2]/rand_old[0,1])**2 + (rand_old[:,2]/rand_old[:,1])**2)
rand_new_speedup_std = rand_new_speedup * \
  np.sqrt((rand_new[0,2]/rand_new[0,1])**2 + (rand_new[:,2]/rand_new[:,1])**2)
reg_old_speedup_std = reg_old_speedup * \
  np.sqrt((reg_old[0,2]/reg_old[0,1])**2 + (reg_old[:,2]/reg_old[:,1])**2)
reg_new_speedup_std = reg_new_speedup * \
  np.sqrt((reg_new[0,2]/reg_new[0,1])**2 + (reg_new[:,2]/reg_new[:,1])**2)

num_threads = len(rand_old)

fig, ax = pl.subplots(2, 2, sharex = True, sharey = "row")

ax[0][0].errorbar(rand_old[:,0], rand_old[:,1], yerr = rand_old[:,2], fmt = '.',
                  label = "old")
ax[0][0].errorbar(rand_new[:,0], rand_new[:,1], yerr = rand_new[:,2], fmt = '.',
                  label = "new")
ax[0][0].set_xscale("log", basex = 2)
ax[0][0].set_title("random")
ax[0][0].legend(loc = "best")
ax[0][0].set_ylabel("total time (s)")
ax[0][0].set_xlim(0.9, num_threads*1.1)

ax[0][1].errorbar(reg_old[:,0], reg_old[:,1], yerr = reg_old[:,2], fmt = '.',
                  label = "old")
ax[0][1].errorbar(reg_new[:,0], reg_new[:,1], yerr = reg_new[:,2], fmt = '.',
                  label = "new")
ax[0][1].set_xscale("log", basex = 2)
ax[0][1].set_title("regular")
ax[0][1].legend(loc = "best")

ax[1][0].plot([1, num_threads], [1, num_threads], "k-")
ax[1][0].errorbar(rand_old[:,0], rand_old_speedup, yerr = rand_old_speedup_std,
                  fmt = '.', label = "old")
ax[1][0].errorbar(rand_new[:,0], rand_new_speedup, yerr = rand_new_speedup_std,
                  fmt = '.', label = "new")
ax[1][0].set_xscale("log", basex = 2)
ax[1][0].set_yscale("log", basey = 2)
ax[1][0].legend(loc = "best")
ax[1][0].set_xlabel("number of threads")
ax[1][0].set_ylabel("speedup")
ax[1][0].set_ylim(0.9, num_threads*1.1)

ax[1][1].plot([1, num_threads], [1, num_threads], "k-")
ax[1][1].errorbar(reg_old[:,0], reg_old_speedup, yerr = reg_old_speedup_std, 
                  fmt = '.', label = "old")
ax[1][1].errorbar(reg_new[:,0], reg_new_speedup, yerr = reg_new_speedup_std,
                  fmt = '.', label = "new")
ax[1][1].set_xscale("log", basex = 2)
ax[1][1].set_yscale("log", basey = 2)
ax[1][1].legend(loc = "best")
ax[1][1].set_xlabel("number of threads")

pl.tight_layout()
pl.savefig("timeVoronoiGrids_scaling.png")
