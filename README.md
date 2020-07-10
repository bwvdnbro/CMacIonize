![Automated build and unit tests](https://github.com/bwvdnbro/CMacIonize/workflows/Automated%20build%20and%20unit%20tests/badge.svg)

# CMacIonize
C++ Monte Carlo photoionization and radiation hydrodynamics code.

 - [Obtaining and compiling the code](#compilation)
 - [Code dependencies](#dependencies)
 - [Documentation](#documentation)
 - [Support](#support)

## Compilation

First, clone the repository using
```
git clone https://github.com/bwvdnbro/CMacIonize.git
```
This creates a new folder named `CMacIonize`. Enter this folder and create a new
`build` folder:
```
cd CMacIonize
mkdir build
```
Enter the `build` folder and call CMake:
```
cd build
cmake ..
```
Now compile the program using `make`.

Documentation can be generated using `make doc`.

Unit tests can be run using `make check`.

A number of benchmark problems is provided in the folder
`build/rundir/benchmark`. Each benchmark problem consists of at least a
parameter file that can be used to run the benchmark, and a `.txt` file
containing some more information about the benchmark problem.

## Dependencies

The only thing you need to compile and run CMacIonize is the CMake build system
and a C++ compiler. However, some plugins for reading and writing files require
the HDF5 library as well. The CMake build system automatically tries to locate
a working HDF5 installation on your system. If this fails, the program will
still compile, but the following functionality will not work:
  - reading a density field from a Gadget snapshot file
  - reading a density field from a FLASH snapshot file
  - reading ionizing sources from a Gadget snapshot file
  - writing Gadget output files. Since this is currently the only supported
    output format, this means output will not work, and the program is in fact
    useless.

To install the HDF5 (development) libraries on a Linux system, use
```
sudo apt-get install libhdf5-dev
```

On Mac OS X, the easiest way to install the HDF5 (development) libraries is a
manual installation. Get the latest HDF5 release tarball from the
[HDF5 website](https://support.hdfgroup.org/HDF5/release/obtainsrc.html)
(make sure to choose a distribution with UNIX line endings). Extract the
tarball in a folder of choice, and configure and build the program:
```
./configure
make
make install
```
Depending on where the HDF5 library is installed, you might need to tell CMake
where to find HDF5: locate the folder containing `libhdf5.so` (or similar).
This folder should be named `lib`, and should be part of a parent folder that
also contains folders named `include` and `bin`. You should set the environment
variable `HDF5_ROOT` to the name of this folder:
```
export HDF5_ROOT=/name/of/the/folder
```
Check the messages generated when you run CMake, as it will tell you whether or
not a working HDF5 development installation was found.

CMacIonize uses Doxygen to generate code documentation, but this is in no way
necessary to compile or run the program. CMake will automatically try to find
a Doxygen installation. If it finds one, you can generate documentation using
`make doc`. If not, `make doc` will not work.

Some scripts used for analyzing benchmark results, or for writing files used
by the unit tests, require Python, and the `numpy`, `matplotlib` and `h5py`
libraries. However, all files that are necessary for the unit tests are also
included in the repository (even if the scripts to generate them are), so
it is still possible to compile and run the code without Python.

## Documentation

The CMacIonize code contains a full inline documentation using Doxygen. A recent version of this documentation is available from [an online mirror](https://users.ugent.be/~bwvdnbro/CMacIonize/). A small number of [online tutorials](https://bwvdnbro.github.io/CMacIonize/tutorials/) is available from [the CMacIonize webpage](https://bwvdnbro.github.io/CMacIonize/). Other sources of documentation include:
 - the journal article describing the original photoionization code on which CMacIonize was based, [Wood, Mathis & Ercolano (2004)](https://ui.adsabs.harvard.edu/abs/2004MNRAS.348.1337W/abstract)
 - the journal article describing CMacIonize 1.0, [Vandenbroucke & Wood (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26C....23...40V/abstract)
 - the journal article describing the task-based algorithm underlying CMacIonize 2.0, [Vandenbroucke & Camps (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200615147V/abstract)

## Support

If you have any issues obtaining, compiling or running the code, please consult the [wiki](https://github.com/bwvdnbro/CMacIonize/wiki) or the [online tutorials](https://bwvdnbro.github.io/CMacIonize/tutorials/), or let us know by email (bert.vandenbroucke@gmail.com). A [Slack workspace](https://cmacionize.slack.com) is available for regular users and contributors upon invitation. Try to use github issue were applicable to
 - [submit bug reports](https://github.com/bwvdnbro/CMacIonize/issues/new?template=bug.md)
 - [request new code features](https://github.com/bwvdnbro/CMacIonize/issues/new?template=enhancement.md)
