import numpy as np
import h5py
import re

# precompiled regular expressions for reading in int and length arrays
intarray = re.compile("\[([0-9]*), ([0-9]*), ([0-9]*)\]")
lengtharray = re.compile(
    "\[([\-0-9\.e\+]*) m, ([\-0-9\.e\+]*) m, ([\-0-9\.e\+]*) m\]"
)

##
# @brief Read quantities from a distributed snapshot. This method takes an
# optional additional conversion function that can be used to convert the
# quantities.
#
# The quantities are read in chunks, with each chunk corresponding to a single
# subgrid in the distributed grid. If @f$N@f$ is the number of quantities and
# @f$s@f$ is the size of a chunk, then the datablock after reading has shape
# @f$(N, s)@f$. The optional conversion function takes the @f$N@f$ input
# quantities and converts them to @f$N'@f$ output quantities:
# @f[
#    f: (N, s) \rightarrow{} (N', s)
# @f]
# If no conversion function is provided, @f$N' = N@f$.
#
# The output is a data cube with shape @f$(N', N_x, N_y, N_z)@f$, where
# @f$N_x@f$, @f$N_y@f$ and @f$N_z@f$ are the number of cells in the grid in each
# coordinate direction. If @f$N'=1@f$, the first dimension of the image cube is
# squeezed and a 3D array is returned.
#
# @param filename Name of the snapshot file.
# @param quantities str or list of str containing the names of all quantities
# to read.
# @param conversion_function (optional) Conversion function applied to the
# quantities after reading. Default: no conversion, quantities are returned as
# they are read.
# @param output_coordinates (optional) Output the cell coordinates? Default is
# false.
# @return 1 or 4 arrays: coordinates of the cells (3D, only if
# output_coordinates is true), data cube containing the quantities (3D or 4D,
# depending on the number of quantities or the number of quantities returned by
# the conversion function).
##
def read_quantities(
    filename,
    quantities,
    conversion_function=lambda x: x,
    output_coordinates=False,
):

    # convert quantities to a list if it is not one already
    if not isinstance(quantities, list):
        quantities = [quantities]

    # open the snapshot file
    file = h5py.File(filename, "r")

    # make sure we are dealing with a distributed grid
    if (
        not "DensitySubGridCreator:number of subgrids"
        in file["/Parameters"].attrs
    ):
        print("Error: not a distributed snapshot!")
        exit(1)

    # read grid parameters
    res_x, res_y, res_z = intarray.search(
        file["/Parameters"].attrs["DensityGrid:number of cells"].decode("utf-8")
    ).groups()
    res = np.array([int(res_x), int(res_y), int(res_z)], dtype=np.uint32)
    ch_x, ch_y, ch_z = intarray.search(
        file["/Parameters"]
        .attrs["DensitySubGridCreator:number of subgrids"]
        .decode("utf-8")
    ).groups()
    chunks = np.array([int(ch_x), int(ch_y), int(ch_z)], dtype=np.uint32)
    chunksize = res // chunks

    if output_coordinates:
        # read box dimensions
        ax, ay, az = lengtharray.search(
            file["/Parameters"].attrs["SimulationBox:anchor"].decode("utf-8")
        ).groups()
        sx, sy, sz = lengtharray.search(
            file["/Parameters"].attrs["SimulationBox:sides"].decode("utf-8")
        ).groups()
        box_anchor = np.array(
            [float(ax), float(ay), float(az)], dtype=np.float32
        )
        box_sides = np.array(
            [float(sx), float(sy), float(sz)], dtype=np.float32
        )

    # set up an empty data cube (we will initialise it once the output shape of
    # conversion_function() is known
    qgrid = None

    # control variables for the reading loop
    startchunk = 0
    endchunk = chunksize[0] * chunksize[1] * chunksize[2]
    chunk_length = chunksize[0] * chunksize[1] * chunksize[2]
    max_chunk = res[0] * res[1] * res[2]
    chunk_shape = None
    ix = 0
    iy = 0
    iz = 0
    # reading loop
    while endchunk <= max_chunk:
        # create a data cube for the chunk data
        qs = np.zeros((len(quantities), endchunk - startchunk))
        # read the requested quantities
        for iq in range(len(quantities)):
            qs[iq] = file["/PartType0/" + quantities[iq]][startchunk:endchunk]

        # apply the conversion function
        qs = conversion_function(qs)

        # qgrid.shape[0] now contains the number of output quantities
        # use this to initialise the data cube and to determine the shape of
        # each converted chunk
        if qgrid is None:
            if len(qs.shape) == 1:
                qgrid = np.zeros((1, res[0], res[1], res[2]))
                chunk_shape = (1, chunksize[0], chunksize[1], chunksize[2])
            else:
                qgrid = np.zeros((qs.shape[0], res[0], res[1], res[2]))
                chunk_shape = (
                    qs.shape[0],
                    chunksize[0],
                    chunksize[1],
                    chunksize[2],
                )

        # put the values in the appropriate cells
        qgrid[
            :,
            ix : ix + chunksize[0],
            iy : iy + chunksize[1],
            iz : iz + chunksize[2],
        ] = qs.reshape(chunk_shape)

        # update the read loop control variables
        startchunk += chunk_length
        endchunk += chunk_length
        iz += chunksize[2]
        if iz == res[2]:
            iz = 0
            iy += chunksize[1]
            if iy == res[1]:
                iy = 0
                ix += chunksize[0]

    if output_coordinates:
        # generate the cell coordinate arrays:
        # first generate 1D arrays for the boundaries of the pixels
        ixgrid = np.linspace(
            box_anchor[0], box_anchor[0] + box_sides[0], res[0] + 1
        )
        iygrid = np.linspace(
            box_anchor[1], box_anchor[1] + box_sides[1], res[1] + 1
        )
        izgrid = np.linspace(
            box_anchor[2], box_anchor[2] + box_sides[2], res[2] + 1
        )
        # now compute the center of each pixel
        ixgrid = 0.5 * (ixgrid[1:] + ixgrid[:-1])
        iygrid = 0.5 * (iygrid[1:] + iygrid[:-1])
        izgrid = 0.5 * (izgrid[1:] + izgrid[:-1])
        # set up a 3D grid
        xgrid, ygrid, zgrid = np.meshgrid(ixgrid, iygrid, izgrid, indexing="ij")

        # return and exit the function
        return xgrid, ygrid, zgrid, qgrid.squeeze()
    else:
        # return and exit the function
        return qgrid.squeeze()


##
# @brief Read quantities from a distributed snapshot and integrate them out
# along the given coordinate axis. This method takes an optional additional
# conversion function that can be used to convert the quantities before the
# integration.
#
# The quantities are read in chunks, with each chunk corresponding to a single
# subgrid in the distributed grid. If @f$N@f$ is the number of quantities and
# @f$s@f$ is the size of a chunk, then the datablock after reading has shape
# @f$(N, s)@f$. The optional conversion function takes the @f$N@f$ input
# quantities and converts them to @f$N'@f$ output quantities:
# @f[
#    f: (N, s) \rightarrow{} (N', s)
# @f]
# If no conversion function is provided, @f$N' = N@f$.
#
# The output is an image data cube with shape @f$(N', N_i, N_j)@f$, where
# @f$N_i@f$ and @f$N_j@f$ are the number of cells in the grid in the plane
# perpendicular to the integration axis. If @f$N'=1@f$, the first dimension of
# the image cube is squeezed and a 2D array is returned.
#
# @param filename Name of the snapshot file.
# @param quantities str or list of str containing the names of all quantities
# to read.
# @param conversion_function (optional) Conversion function applied to the
# quantities after reading but before integrating along the line of sight.
# Default: no conversion, quantities are integrated as they are read.
# @param axis (optional) Axis along which to integrate. Defaults to the z-axis
# (axis 2).
# @param output_coordinates (optional) Output the pixel coordinates? Default is
# false.
# @return 1 or 3 arrays: coordinates of the image pixels (2D, only if
# output_coordinates is true), image cube containing the integrated quantities
# (2D or 3D, depending on the number of quantities or the number of quantities
# returned by the conversion function).
##
def read_integrated_quantities(
    filename,
    quantities,
    conversion_function=lambda x: x,
    axis=2,
    output_coordinates=False,
):

    # convert quantities to a list if it is not one already
    if not isinstance(quantities, list):
        quantities = [quantities]

    # open the snapshot file
    file = h5py.File(filename, "r")

    # make sure we are dealing with a distributed grid
    if (
        not "DensitySubGridCreator:number of subgrids"
        in file["/Parameters"].attrs
    ):
        print("Error: not a distributed snapshot!")
        exit(1)

    # read grid parameters
    res_x, res_y, res_z = intarray.search(
        file["/Parameters"].attrs["DensityGrid:number of cells"].decode("utf-8")
    ).groups()
    res = np.array([int(res_x), int(res_y), int(res_z)], dtype=np.uint32)
    ch_x, ch_y, ch_z = intarray.search(
        file["/Parameters"]
        .attrs["DensitySubGridCreator:number of subgrids"]
        .decode("utf-8")
    ).groups()
    chunks = np.array([int(ch_x), int(ch_y), int(ch_z)], dtype=np.uint32)
    chunksize = res // chunks

    if output_coordinates:
        # read box dimensions
        ax, ay, az = lengtharray.search(
            file["/Parameters"].attrs["SimulationBox:anchor"].decode("utf-8")
        ).groups()
        sx, sy, sz = lengtharray.search(
            file["/Parameters"].attrs["SimulationBox:sides"].decode("utf-8")
        ).groups()
        box_anchor = np.array(
            [float(ax), float(ay), float(az)], dtype=np.float32
        )
        box_sides = np.array(
            [float(sx), float(sy), float(sz)], dtype=np.float32
        )

    # determine the dimensions of the image
    if axis == 0:
        image_shape = (res[1], res[2])
    elif axis == 1:
        image_shape = (res[0], res[2])
    else:
        image_shape = (res[0], res[1])

    # set up an empty image data cube (we will initialise it once the output
    # shape of conversion_function() is known
    qgrid = None

    # control variables for the reading loop
    startchunk = 0
    endchunk = chunksize[0] * chunksize[1] * chunksize[2]
    chunk_length = chunksize[0] * chunksize[1] * chunksize[2]
    max_chunk = res[0] * res[1] * res[2]
    chunk_shape = None
    ix = 0
    iy = 0
    iz = 0
    # reading loop
    while endchunk <= max_chunk:
        # create a data cube for the chunk data
        qs = np.zeros((len(quantities), endchunk - startchunk))
        # read the requested quantities
        for iq in range(len(quantities)):
            qs[iq] = file["/PartType0/" + quantities[iq]][startchunk:endchunk]

        # apply the conversion function
        qs = conversion_function(qs)

        # qgrid.shape[0] now contains the number of output quantities
        # use this to initialise the data cube and to determine the shape of
        # each converted chunk
        if qgrid is None:
            if len(qs.shape) == 1:
                qgrid = np.zeros((1, image_shape[0], image_shape[1]))
                chunk_shape = (1, chunksize[0], chunksize[1], chunksize[2])
            else:
                qgrid = np.zeros((qs.shape[0], image_shape[0], image_shape[1]))
                chunk_shape = (
                    qs.shape[0],
                    chunksize[0],
                    chunksize[1],
                    chunksize[2],
                )

        # perform the integration and add the result to the appropriate pixels
        if axis == 0:
            qgrid[
                :, iy : iy + chunksize[1], iz : iz + chunksize[2]
            ] += qs.reshape(chunk_shape).sum(axis=axis + 1)
        elif axis == 1:
            qgrid[
                :, ix : ix + chunksize[0], iz : iz + chunksize[2]
            ] += qs.reshape(chunk_shape).sum(axis=axis + 1)
        else:
            qgrid[
                :, ix : ix + chunksize[0], iy : iy + chunksize[1]
            ] += qs.reshape(chunk_shape).sum(axis=axis + 1)

        # update the read loop control variables
        startchunk += chunk_length
        endchunk += chunk_length
        iz += chunksize[2]
        if iz == res[2]:
            iz = 0
            iy += chunksize[1]
            if iy == res[1]:
                iy = 0
                ix += chunksize[0]

    # close the file, we are done with it
    file.close()

    if output_coordinates:
        # generate the pixel coordinate arrays:
        # first generate 1D arrays for the boundaries of the pixels
        if axis == 0:
            ixgrid = np.linspace(
                box_anchor[1], box_anchor[1] + box_sides[1], res[1] + 1
            )
            iygrid = np.linspace(
                box_anchor[2], box_anchor[2] + box_sides[2], res[2] + 1
            )
        elif axis == 1:
            ixgrid = np.linspace(
                box_anchor[0], box_anchor[0] + box_sides[0], res[0] + 1
            )
            iygrid = np.linspace(
                box_anchor[2], box_anchor[2] + box_sides[2], res[2] + 1
            )
        else:
            ixgrid = np.linspace(
                box_anchor[0], box_anchor[0] + box_sides[0], res[0] + 1
            )
            iygrid = np.linspace(
                box_anchor[1], box_anchor[1] + box_sides[1], res[1] + 1
            )
        # now compute the center of each pixel
        ixgrid = 0.5 * (ixgrid[1:] + ixgrid[:-1])
        iygrid = 0.5 * (iygrid[1:] + iygrid[:-1])
        # set up a 2D grid
        xgrid, ygrid = np.meshgrid(ixgrid, iygrid, indexing="ij")

        # return and exit the function
        return xgrid, ygrid, qgrid.squeeze()
    else:
        # return and exit the function
        return qgrid.squeeze()


##
# @brief Read quantities from a distributed snapshot and slice them
# perpendicular to the given axis at the given intercept. This method takes an
# optional additional conversion function that can be used to convert the
# quantities before the slicing.
#
# The quantities are read in chunks, with each chunk corresponding to a single
# subgrid in the distributed grid. If @f$N@f$ is the number of quantities and
# @f$s@f$ is the size of a chunk, then the datablock after reading has shape
# @f$(N, s)@f$. The optional conversion function takes the @f$N@f$ input
# quantities and converts them to @f$N'@f$ output quantities:
# @f[
#    f: (N, s) \rightarrow{} (N', s)
# @f]
# If no conversion function is provided, @f$N' = N@f$.
#
# The output is an image data cube with shape @f$(N', N_i, N_j)@f$, where
# @f$N_i@f$ and @f$N_j@f$ are the number of cells in the grid in the plane
# perpendicular to the slice axis. If @f$N'=1@f$, the first dimension of the
# image cube is squeezed and a 2D array is returned.
#
# @param filename Name of the snapshot file.
# @param quantities str or list of str containing the names of all quantities
# to read.
# @param conversion_function (optional) Conversion function applied to the
# quantities after reading but before slicing perpendicular to the line of
# sight. Default: no conversion, quantities are sliced as they are read.
# @param axis (optional) Axis perpendicular to which to integrate. Defaults to
# the z-axis (axis 2).
# @param intercept (optional) Intercept along the slice axis, in SI units.
# Default is the origin (z = 0).
# @param output_coordinates (optional) Output the pixel coordinates? Default is
# false.
# @return 1 or 3 arrays: coordinates of the image pixels (2D, only if
# output_coordinates is true), image cube containing the sliced quantities
# (2D or 3D, depending on the number of quantities or the number of quantities
# returned by the conversion function).
##
def read_sliced_quantities(
    filename,
    quantities,
    conversion_function=lambda x: x,
    axis=2,
    intercept=0,
    output_coordinates=False,
):

    # convert quantities to a list if it is not one already
    if not isinstance(quantities, list):
        quantities = [quantities]

    # open the snapshot file
    file = h5py.File(filename, "r")

    # make sure we are dealing with a distributed grid
    if (
        not "DensitySubGridCreator:number of subgrids"
        in file["/Parameters"].attrs
    ):
        print("Error: not a distributed snapshot!")
        exit(1)

    # read grid parameters
    res_x, res_y, res_z = intarray.search(
        file["/Parameters"].attrs["DensityGrid:number of cells"].decode("utf-8")
    ).groups()
    res = np.array([int(res_x), int(res_y), int(res_z)], dtype=np.uint32)
    ch_x, ch_y, ch_z = intarray.search(
        file["/Parameters"]
        .attrs["DensitySubGridCreator:number of subgrids"]
        .decode("utf-8")
    ).groups()
    chunks = np.array([int(ch_x), int(ch_y), int(ch_z)], dtype=np.uint32)
    chunksize = res // chunks

    # read box dimensions
    ax, ay, az = lengtharray.search(
        file["/Parameters"].attrs["SimulationBox:anchor"].decode("utf-8")
    ).groups()
    sx, sy, sz = lengtharray.search(
        file["/Parameters"].attrs["SimulationBox:sides"].decode("utf-8")
    ).groups()
    box_anchor = np.array([float(ax), float(ay), float(az)], dtype=np.float32)
    box_sides = np.array([float(sx), float(sy), float(sz)], dtype=np.float32)

    # check that the intercept is within the box
    if (
        intercept < box_anchor[axis]
        or intercept > box_anchor[axis] + box_sides[axis]
    ):
        print("Error: intercept value outside box range!")
        exit(1)
    # determine the index of the slice within the cell grid
    slice_index = int(
        np.floor((intercept - box_anchor[axis]) * res[axis] / box_sides[axis])
    )

    # determine the dimensions of the image
    if axis == 0:
        image_shape = (res[1], res[2])
    elif axis == 1:
        image_shape = (res[0], res[2])
    else:
        image_shape = (res[0], res[1])

    # set up an empty image data cube (we will initialise it once the output
    # shape of conversion_function() is known
    qgrid = None

    # control variables for the reading loop
    startchunk = 0
    endchunk = chunksize[0] * chunksize[1] * chunksize[2]
    chunk_length = chunksize[0] * chunksize[1] * chunksize[2]
    max_chunk = res[0] * res[1] * res[2]
    chunk_shape = None
    ix = 0
    iy = 0
    iz = 0
    # reading loop
    while endchunk <= max_chunk:

        # check if the current chunk overlaps with the slice
        if axis == 0:
            in_slice = (ix <= slice_index) and (slice_index < ix + chunksize[0])
        elif axis == 1:
            in_slice = (iy <= slice_index) and (slice_index < iy + chunksize[1])
        else:
            in_slice = (iz <= slice_index) and (slice_index < iz + chunksize[2])

        if in_slice:
            # create a data cube for the chunk data
            qs = np.zeros((len(quantities), endchunk - startchunk))
            # read the requested quantities
            for iq in range(len(quantities)):
                qs[iq] = file["/PartType0/" + quantities[iq]][
                    startchunk:endchunk
                ]

            # apply the conversion function
            qs = conversion_function(qs)

            # qgrid.shape[0] now contains the number of output quantities
            # use this to initialise the data cube and to determine the shape of
            # each converted chunk
            if qgrid is None:
                if len(qs.shape) == 1:
                    qgrid = np.zeros((1, image_shape[0], image_shape[1]))
                    chunk_shape = (1, chunksize[0], chunksize[1], chunksize[2])
                else:
                    qgrid = np.zeros(
                        (qs.shape[0], image_shape[0], image_shape[1])
                    )
                    chunk_shape = (
                        qs.shape[0],
                        chunksize[0],
                        chunksize[1],
                        chunksize[2],
                    )

            # slice the chunk and put the results in the right image pixel
            if axis == 0:
                qgrid[
                    :, iy : iy + chunksize[1], iz : iz + chunksize[2]
                ] = qs.reshape(chunk_shape)[:, slice_index - ix, :, :]
            elif axis == 1:
                qgrid[
                    :, ix : ix + chunksize[0], iz : iz + chunksize[2]
                ] = qs.reshape(chunk_shape)[:, :, slice_index - iy, :]
            else:
                qgrid[
                    :, ix : ix + chunksize[0], iy : iy + chunksize[1]
                ] = qs.reshape(chunk_shape)[:, :, :, slice_index - iz]

        # update the read loop control variables
        startchunk += chunk_length
        endchunk += chunk_length
        iz += chunksize[2]
        if iz == res[2]:
            iz = 0
            iy += chunksize[1]
            if iy == res[1]:
                iy = 0
                ix += chunksize[0]

    # close the file, we are done with it
    file.close()

    if output_coordinates:
        # generate the pixel coordinate arrays:
        # first generate 1D arrays for the boundaries of the pixels
        if axis == 0:
            ixgrid = np.linspace(
                box_anchor[1], box_anchor[1] + box_sides[1], res[1] + 1
            )
            iygrid = np.linspace(
                box_anchor[2], box_anchor[2] + box_sides[2], res[2] + 1
            )
        elif axis == 1:
            ixgrid = np.linspace(
                box_anchor[0], box_anchor[0] + box_sides[0], res[0] + 1
            )
            iygrid = np.linspace(
                box_anchor[2], box_anchor[2] + box_sides[2], res[2] + 1
            )
        else:
            ixgrid = np.linspace(
                box_anchor[0], box_anchor[0] + box_sides[0], res[0] + 1
            )
            iygrid = np.linspace(
                box_anchor[1], box_anchor[1] + box_sides[1], res[1] + 1
            )
        # now compute the center of each pixel
        ixgrid = 0.5 * (ixgrid[1:] + ixgrid[:-1])
        iygrid = 0.5 * (iygrid[1:] + iygrid[:-1])
        # set up a 2D grid
        xgrid, ygrid = np.meshgrid(ixgrid, iygrid, indexing="ij")

        # return and exit the function
        return xgrid, ygrid, qgrid.squeeze()
    else:
        # return and exit the function
        return qgrid.squeeze()


##
# @brief Unit test.
#
# Takes the name of a snaphot file as input and generates an image of the
# density of neutral density using the read, integrate or slice method. Further
# arguments allow to test the method versions with and without coordinate
# arrays.
##
if __name__ == "__main__":
    import argparse
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pl

    pl.rcParams["text.usetex"] = True

    # command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--file", "-f", action="store", required=True)
    argparser.add_argument("--output", "-o", action="store", required=True)
    argparser.add_argument(
        "--method",
        "-m",
        action="store",
        choices=["read", "integrate", "slice"],
        default="slice",
    )
    argparser.add_argument("--neutral", "-n", action="store_true")
    argparser.add_argument("--coordinates", "-c", action="store_true")
    args = argparser.parse_args()

    # function used to compute the neutral density based on the density and
    # the neutral fraction
    def neutral_density_function(quantities):
        return quantities[0] * quantities[1]

    # call the test function with the appropriate arguments
    if args.method == "read":
        if args.neutral:
            result = read_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                neutral_density_function,
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, zgrid, rhogrid = result
                xgrid = xgrid[:, :, 0]
                ygrid = ygrid[:, :, 0]
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 3])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
        else:
            result = read_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, zgrid, rhogrid = result
                xgrid = xgrid[:, :, 0]
                ygrid = ygrid[:, :, 0]
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 3])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
            rhogrid = rhogrid[0]
        rhogrid = rhogrid.sum(axis=2)
    elif args.method == "integrate":
        if args.neutral:
            result = read_integrated_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                neutral_density_function,
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, rhogrid = result
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 1])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
        else:
            result = read_integrated_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, rhogrid = result
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 1])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
            rhogrid = rhogrid[0]
    elif args.method == "slice":
        if args.neutral:
            result = read_sliced_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                neutral_density_function,
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, rhogrid = result
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 1])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
        else:
            result = read_sliced_quantities(
                args.file,
                ["NumberDensity", "NeutralFractionH"],
                output_coordinates=args.coordinates,
            )
            if args.coordinates:
                xgrid, ygrid, rhogrid = result
            else:
                rhogrid = result
                xgrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 2])
                ygrid = np.arange(0, rhogrid.shape[len(rhogrid.shape) - 1])
                xgrid, ygrid = np.meshgrid(xgrid, ygrid, indexing="ij")
            rhogrid = rhogrid[0]

    # plot the result
    pl.contourf(xgrid, ygrid, np.log10(rhogrid), 500)
    pl.savefig(args.output, dpi=300, bbox_inches="tight")
