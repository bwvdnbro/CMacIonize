#! /usr/bin/python

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
# @file bondi.py
#
# @brief Read the final Bondi snapshot file and plot it together with the
# analytic solution.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# load some modules
import numpy as np
import h5py

# overwrite the default matplotlib backend to support plotting without an active
# graphics environment
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

# load the lambert W function, as the analytic solution uses it
import scipy.special.lambertw as lambertw

# load statistical tools to bin simulation data
import scipy.stats as stats

# overwrite the default matplotlib settings to get a nicer figure
pl.rcParams["figure.figsize"] = (8, 10)
pl.rcParams["text.usetex"] = True

##
# @brief Bin the given quantity in the given number of equally spaced bins.
#
# @param x Quantity positions.
# @param q Quantity to bin.
# @param nbin Number of spatial bins.
# @return Bin positions, mean quantity in each bin, standard deviation within
# each bin.
##
def bin_quantity(x, q, nbin):
    # we let binned_statistic figure out the bins and then reuse them for the
    # standard deviation mean
    qbin, bins, _ = stats.binned_statistic(x, q, statistic="mean", bins=nbin)
    q2bin, _, _ = stats.binned_statistic(x, q ** 2, statistic="mean", bins=bins)
    qsigma = np.sqrt(q2bin - qbin ** 2)
    return 0.5 * (bins[1:] + bins[:-1]), qbin, qsigma


# we will plot distances in AU rather than m
au_in_si = 1.495978707e11  # m

## Bondi parameters
# input unit parameters
G_in_si = 6.67408e-11
k_in_si = 1.38064852e-23
mH_in_si = 1.674e-27
solar_mass_in_si = 1.9891e30

# input parameters
# physical
mass_point_mass = 18.0 * solar_mass_in_si
T_n = 500
bondi_rho_n = 1.0e-16

# derived parameters
cs2_n = T_n * k_in_si / mH_in_si
bondi_r_n = 0.5 * G_in_si * mass_point_mass / cs2_n
cs_n = np.sqrt(cs2_n)

## compute the analytic solution

##
# @brief Get the analytic Bondi solution for the given positions.
#
# @param r Positions.
# @return Density, velocity and pressure.
##
def bondi_analytic(r):
    u = bondi_r_n / r
    omega = -u ** 4 * np.exp(3.0 - 4.0 * u)
    v_a = np.where(
        r < bondi_r_n,
        -cs_n * np.sqrt(-lambertw(omega, -1).real),
        -cs_n * np.sqrt(-lambertw(omega, 0).real),
    )
    rho_a = -bondi_rho_n * bondi_r_n ** 2 * cs_n / (r ** 2 * v_a)
    P_a = cs2_n * rho_a

    return rho_a, v_a, P_a


r_a = np.arange(
    10.1 * au_in_si, 50.0 * np.sqrt(3.0) * au_in_si, 0.01 * au_in_si
)
rho_a, v_a, P_a = bondi_analytic(r_a)

# open the final snapshot filoe
file = h5py.File("bondi_001.hdf5", "r")

# read the box size, coordinates and hydro quantities
box = file["/Header"].attrs["BoxSize"]
coords = file["/PartType0/Coordinates"]
rho = file["/PartType0/Density"]
vs = file["/PartType0/Velocities"]
P = file["/PartType0/Pressure"]

# convert coordinates to radii
radius = np.sqrt(
    (coords[:, 0] - 0.5 * box[0]) ** 2
    + (coords[:, 1] - 0.5 * box[1]) ** 2
    + (coords[:, 2] - 0.5 * box[2]) ** 2
)
# convert velocities to radial velocity
v = (
    vs[:, 0] * (coords[:, 0] - 0.5 * box[0])
    + vs[:, 1] * (coords[:, 1] - 0.5 * box[1])
    + vs[:, 2] * (coords[:, 2] - 0.5 * box[2])
) / radius

# create the plot
fig, ax = pl.subplots(3, 2, sharex=True)

# plot the quantities:
# density
ax[0][0].semilogy(
    radius / au_in_si, rho[:] * 0.001, "k.", markersize=0.5, alpha=0.1
)
rbin, rhobin, rhosigma = bin_quantity(radius, rho[:], 100)
ax[0][0].errorbar(
    rbin / au_in_si,
    rhobin * 0.001,
    yerr=rhosigma * 0.001,
    label="\\sc{}CMacIonize",
    fmt=".",
)
rind = rbin > 10.5 * au_in_si
rhoref, _, _ = bondi_analytic(rbin[rind])
ax[0][1].plot(
    rbin[rind] / au_in_si,
    (rhobin[rind] - rhoref) / (rhobin[rind] + rhoref),
    ".",
)

# velocity
ax[1][0].plot(radius / au_in_si, v * 0.001, "k.", markersize=0.5, alpha=0.1)
rbin, vbin, vsigma = bin_quantity(radius, v, 100)
ax[1][0].errorbar(rbin / au_in_si, vbin * 0.001, yerr=vsigma * 0.001, fmt=".")
rind = rbin > 10.5 * au_in_si
_, vref, _ = bondi_analytic(rbin[rind])
ax[1][1].plot(
    rbin[rind] / au_in_si, (vbin[rind] - vref) / (vbin[rind] + vref), "."
)

# pressure
ax[2][0].semilogy(radius / au_in_si, P, "k.", markersize=0.5, alpha=0.1)
rbin, Pbin, Psigma = bin_quantity(radius, P[:], 100)
ax[2][0].errorbar(rbin / au_in_si, Pbin, yerr=Psigma, fmt=".")
rind = rbin > 10.5 * au_in_si
_, _, Pref = bondi_analytic(rbin[rind])
ax[2][1].plot(
    rbin[rind] / au_in_si, (Pbin[rind] - Pref) / (Pbin[rind] + Pref), "."
)

# analytic solution
ax[0][0].semilogy(r_a / au_in_si, rho_a * 0.001, label="analytic")
ax[1][0].plot(r_a / au_in_si, v_a * 0.001)
ax[2][0].semilogy(r_a / au_in_si, P_a)

# set axis labels
ax[0][0].set_ylabel("$\\rho{}$ (g cm$^{-3}$)")
ax[1][0].set_ylabel("$v_r$ (km s$^{-1}$)")
ax[2][0].set_xlabel("$r$ (AU)")
ax[2][0].set_ylabel("$P$ (Pa)")

ax[0][1].set_ylabel("$(\\rho{} - \\rho{}_a) / (\\rho{} + \\rho{}_a)$")
ax[1][1].set_ylabel("$(v_r - v_{r,a}) / (v_r + v_{r,a})$")
ax[2][1].set_ylabel("$(P - P_a) / (P + P_a)$")
ax[2][1].set_xlabel("$r$ (AU)")

# add a legend to the first panel
ax[0][0].legend(loc="best")

# add a line to show the inner mask radius
for axrow in ax:
    for a in axrow:
        a.axvline(x=10.0, linestyle="--", color="k", linewidth=0.8)

# finalize and save the plot
pl.tight_layout()
pl.savefig("bondi.png")
