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
# @file gamma_SIV.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Saraph, H. E. & Storey, P. J. 1999, A&AS, 134, 369
# (http://adsabs.harvard.edu/abs/1999A%26AS..134..369S).
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

## load modules
import numpy as np

# for curve_fit
import scipy.optimize as opt

# for plotting (using a backend that does not require a graphics environment)
import matplotlib

matplotlib.use("Agg")
import pylab as pl

# for the fitting curve
from fitting_curve import (
    fitting_curve,
    print_fit_variables,
    initialize_data_values,
    append_data_values,
    print_data_values,
    get_code,
    jacobian_fitting_curve,
    round_parameters,
)

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
    data = {}

    # data from Saraph & Storey (1999), table 5
    T = np.array(
        [
            1000.0,
            1500.0,
            2000.0,
            2500.0,
            3000.0,
            3500.0,
            4000.0,
            5000.0,
            6000.0,
            7000.0,
            8000.0,
            9000.0,
            10000.0,
            11000.0,
            12000.0,
            13000.0,
            14000.0,
            15000.0,
            16000.0,
            18000.0,
            20000.0,
            24000.0,
            28000.0,
            32000.0,
            36000.0,
            40000.0,
        ]
    )
    data = np.array(
        [
            6.89,
            6.93,
            7.11,
            7.32,
            7.52,
            7.70,
            7.85,
            8.10,
            8.27,
            8.39,
            8.47,
            8.52,
            8.55,
            8.56,
            8.55,
            8.54,
            8.51,
            8.48,
            8.44,
            8.36,
            8.28,
            8.09,
            7.92,
            7.76,
            7.61,
            7.47,
        ]
    )

    # initialize the strings for code and value output
    code = ""
    data_values = initialize_data_values()
    # do the curve fitting
    imin = 0
    imax = len(T)
    # fit the curve
    A, _ = opt.curve_fit(
        fitting_curve,
        T[imin:imax],
        data[imin:imax],
        maxfev=1000000,
        method="trf",
        p0=(0.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.0),
        jac=jacobian_fitting_curve,
    )
    A = round_parameters(*A)
    # compute the xi2 difference between the data values (in the fitting
    # interval) and the curve
    xi2 = sum((data[imin:imax] - fitting_curve(T[imin:imax], *A)) ** 2)
    # output some info
    print("Transition: 0 to 1")
    print_fit_variables(*A)
    print("convergence:", xi2)
    print("validity: [", T[imin], ",", T[imax - 1], "]")
    # write the fitting code for this transition
    code += get_code("SIV", "REMOVE_THIS_BLOCK", *A)
    # add the values to the list strings
    append_data_values(data_values, *A)

    # plot the data and fit for visual comparison
    Trange = np.logspace(3.0, 5.0, 100)
    pl.plot(T, data, "k.")
    pl.plot(Trange, fitting_curve(Trange, *A), "r-")
    pl.xlim(0.0, 1.0e5)
    pl.savefig("tmp/SIV_{key}.png".format(key="G0t1"))
    pl.close()

    # save the plot values in separate files
    dfile = open("tmp/SIV_{key}_data.txt".format(key="G0t1"), "w")
    for i in range(len(T)):
        dfile.write("{T}\t{data}\n".format(T=T[i], data=data[i]))
    dfile.close()
    ffile = open("tmp/SIV_{key}_fit.txt".format(key="G0t1"), "w")
    for i in range(len(Trange)):
        ffile.write(
            "{T}\t{fit}\n".format(T=Trange[i], fit=fitting_curve(Trange[i], *A))
        )
    ffile.close()

    # output the code to put into the LineCoolingData constructor
    print("code:")
    print(code)
    # output the values to put in atom4.dat in Kenny's code (to update the
    # reference values for the unit tests)
    print_data_values(data_values)
