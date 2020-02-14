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
# @file gamma_SII.py
#
# @brief Script that fits curves to the velocity-averaged collision strength
# data from Tayal, S. S. & Zatsarinny, O. 2010, ApJS, 188, 32
# (http://adsabs.harvard.edu/abs/2010ApJS..188...32T).
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

# dictionary that links abbreviated transition names to the full names used in
# LineCoolingData
transitions = {
    "G0t1": "TRANSITION_0_to_1",
    "G0t2": "TRANSITION_0_to_2",
    "G0t3": "TRANSITION_0_to_3",
    "G0t4": "TRANSITION_0_to_4",
    "G1t2": "TRANSITION_1_to_2",
    "G1t3": "TRANSITION_1_to_3",
    "G1t4": "TRANSITION_1_to_4",
    "G2t3": "TRANSITION_2_to_3",
    "G2t4": "TRANSITION_2_to_4",
    "G3t4": "TRANSITION_3_to_4",
}

# main function: computes fits to the data and plots the data and fits for
# visual comparison
# the fitted curve coefficients are printed to the stdout
if __name__ == "__main__":
    data = {}

    # data from Tayal & Zatsarinny (2010), table 5
    T = np.array(
        [
            5000.0,
            7000.0,
            10000.0,
            15000.0,
            20000.0,
            25000.0,
            30000.0,
            40000.0,
            50000.0,
            70000.0,
            100000.0,
        ]
    )
    # 4S3/2 to 2D3/2
    data["G0t1"] = np.array(
        [
            2.66e0,
            2.62e0,
            2.56e0,
            2.48e0,
            2.41e0,
            2.35e0,
            2.30e0,
            2.20e0,
            2.10e0,
            1.93e0,
            1.71e0,
        ]
    )
    # 4S3/2 to 2D5/2
    data["G0t2"] = np.array(
        [
            3.98e0,
            3.91e0,
            3.83e0,
            3.71e0,
            3.61e0,
            3.52e0,
            3.44e0,
            3.30e0,
            3.16e0,
            2.90e0,
            2.58e0,
        ]
    )
    # 2D3/2 to 2D5/2
    data["G1t2"] = np.array(
        [
            7.32e0,
            7.13e0,
            6.89e0,
            6.58e0,
            6.35e0,
            6.17e0,
            6.02e0,
            5.76e0,
            5.51e0,
            5.06e0,
            4.48e0,
        ]
    )
    # 2P1/2 to 2P3/2
    data["G3t4"] = np.array(
        [
            1.76e0,
            1.78e0,
            1.80e0,
            1.82e0,
            1.84e0,
            1.85e0,
            1.86e0,
            1.84e0,
            1.80e0,
            1.70e0,
            1.54e0,
        ]
    )
    # 4S3/2 to 2P1/2
    data["G0t3"] = np.array(
        [
            6.86e-1,
            6.94e-1,
            7.04e-1,
            7.17e-1,
            7.27e-1,
            7.33e-1,
            7.36e-1,
            7.30e-1,
            7.15e-1,
            6.70e-1,
            5.97e-1,
        ]
    )
    # 4S3/2 to 2P3/2
    data["G0t4"] = np.array(
        [
            1.38e0,
            1.39e0,
            1.42e0,
            1.44e0,
            1.46e0,
            1.47e0,
            1.48e0,
            1.47e0,
            1.43e0,
            1.34e0,
            1.19e0,
        ]
    )
    # 2D3/2 to 2P1/2
    data["G1t3"] = np.array(
        [
            1.48e0,
            1.48e0,
            1.47e0,
            1.48e0,
            1.48e0,
            1.49e0,
            1.50e0,
            1.50e0,
            1.50e0,
            1.49e0,
            1.45e0,
        ]
    )
    # 2D3/2 to 2P3/2
    data["G1t4"] = np.array(
        [
            2.40e0,
            2.40e0,
            2.39e0,
            2.39e0,
            2.38e0,
            2.38e0,
            2.38e0,
            2.36e0,
            2.33e0,
            2.25e0,
            2.13e0,
        ]
    )
    # 2D5/2 to 2P1/2
    data["G2t3"] = np.array(
        [
            1.79e0,
            1.78e0,
            1.78e0,
            1.77e0,
            1.76e0,
            1.76e0,
            1.75e0,
            1.73e0,
            1.70e0,
            1.63e0,
            1.53e0,
        ]
    )
    # 2D5/2 to 2P3/2
    data["G2t4"] = np.array(
        [
            4.07e0,
            4.07e0,
            4.06e0,
            4.06e0,
            4.07e0,
            4.08e0,
            4.08e0,
            4.08e0,
            4.06e0,
            3.97e0,
            3.83e0,
        ]
    )

    # initialize the strings for code and value output
    code = ""
    data_values = initialize_data_values()
    # do the curve fitting
    for key in sorted(data):
        imin = 0
        imax = len(T)
        # fit the curve
        A, _ = opt.curve_fit(
            fitting_curve,
            T[imin:imax],
            data[key][imin:imax],
            maxfev=100000,
            p0=(0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0),
            jac=jacobian_fitting_curve,
        )
        A = round_parameters(*A)
        # compute the xi2 difference between the data values (in the fitting
        # interval) and the curve
        xi2 = sum((data[key][imin:imax] - fitting_curve(T[imin:imax], *A)) ** 2)
        # output some info
        print("Transition:", key)
        print_fit_variables(*A)
        print("convergence:", xi2)
        print("validity: [", T[imin], ",", T[imax - 1], "]")
        # write the fitting code for this transition
        code += get_code("SII", transitions[key], *A)
        # add the values to the list strings
        append_data_values(data_values, *A)

        # plot the data and fit for visual comparison
        Trange = np.logspace(3.0, 5.0, 100)
        pl.plot(T, data[key], "k.")
        pl.plot(Trange, fitting_curve(Trange, *A), "r-")
        pl.xlim(0.0, 1.0e5)
        pl.savefig("tmp/SII_{key}.png".format(key=key))
        pl.close()

        # save the plot values in separate files
        dfile = open("tmp/SII_{key}_data.txt".format(key=key), "w")
        for i in range(len(T)):
            dfile.write("{T}\t{data}\n".format(T=T[i], data=data[key][i]))
        dfile.close()
        ffile = open("tmp/SII_{key}_fit.txt".format(key=key), "w")
        for i in range(len(Trange)):
            ffile.write(
                "{T}\t{fit}\n".format(
                    T=Trange[i], fit=fitting_curve(Trange[i], *A)
                )
            )
        ffile.close()

    # output the code to put into the LineCoolingData constructor
    print("code:")
    print(code)
    # output the values to put in atom4.dat in Kenny's code (to update the
    # reference values for the unit tests)
    print_data_values(data_values)
