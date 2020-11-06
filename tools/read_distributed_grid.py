import numpy as np
import h5py
import re

# precompiled regular expressions for reading in int and length arrays
intarray = re.compile("\[([0-9]*), ([0-9]*), ([0-9]*)\]")
lengtharray = re.compile(
    "\[([\-0-9\.e\+]*) m, ([\-0-9\.e\+]*) m, ([\-0-9\.e\+]*) m\]"
)


def read_integrated_quantities(
    filename, quantities, conversion_function=lambda x: x, axis=2
):

    file = h5py.File(filename, "r")

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

    if axis == 0:
        image_shape = (len(quantities), res[1], res[2])
    elif axis == 1:
        image_shape = (len(quantities), res[0], res[2])
    else:
        image_shape = (len(quantities), res[0], res[1])

    qgrid = np.zeros(image_shape)

    startchunk = 0
    endchunk = chunksize[0] * chunksize[1] * chunksize[2]
    chunk_length = chunksize[0] * chunksize[1] * chunksize[2]
    max_chunk = res[0] * res[1] * res[2]
    chunk_shape = (len(quantities), chunksize[0], chunksize[1], chunksize[2])
    ix = 0
    iy = 0
    iz = 0
    while endchunk <= max_chunk:
        qs = np.zeros((len(quantities), endchunk - startchunk))
        for iq in range(len(quantities)):
            qs[iq] = file["/PartType0/" + quantities[iq]][startchunk:endchunk]

        qs = conversion_function(qs)

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

        startchunk += chunk_length
        endchunk += chunk_length

        iz += chunksize[2]
        if iz == res[2]:
            iz = 0
            iy += chunksize[1]
            if iy == res[1]:
                iy = 0
                ix += chunksize[0]

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

    ixgrid = 0.5 * (ixgrid[1:] + ixgrid[:-1])
    iygrid = 0.5 * (iygrid[1:] + iygrid[:-1])
    xgrid, ygrid = np.meshgrid(ixgrid, iygrid, indexing="ij")

    return xgrid, ygrid, qgrid.squeeze()


if __name__ == "__main__":
    import argparse
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pl

    pl.rcParams["text.usetex"] = True

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--file", "-f", action="store", required=True)
    argparser.add_argument("--output", "-o", action="store", required=True)
    args = argparser.parse_args()

    def neutral_density_function(quantities):
        output = np.zeros(quantities.shape)
        output[0] = quantities[0] * quantities[1]
        return output

    xgrid, ygrid, rhogrid = read_integrated_quantities(
        args.file,
        ["NumberDensity", "NeutralFractionH"],
        neutral_density_function,
    )

    pl.contourf(xgrid, ygrid, np.log10(rhogrid[0]), 500)

    pl.savefig(args.output, dpi=300, bbox_inches="tight")
