---
layout: tutorial
title: Getting started with CMacIonize
date: 2019-01-11
level: basic
---

{% include toc.html %}

To get started with CMacIonize, you need to do 3 things:

 1. Download the code.
 2. Configure and compile the code.
 3. Run an example.

For this tutorial, we are going to assume that you are working on a 
Unix-based operating system, i.e. either Linux or Mac OSX. CMacIonize 
can also be run on Windows systems, but this is harder to set up and 
might be detailed in a future tutorial.

# Requirements

To obtain and compile a minimal version of CMacIonize, you will need the 
following software and libraries installed on your system:

 * [git](https://git-scm.com/)
 * [CMake](https://cmake.org/)
 * a C/C++ compiler (like [GCC](https://www.gnu.org/software/gcc/) or
   [Clang](https://clang.llvm.org/))

On Linux systems, you should be able to install these using your package 
manager. On Mac OSX, you can install 
[Xcode](https://developer.apple.com/xcode/) to obtain C/C++ compilers 
and you will need to install git and CMake manually. Note that git is 
not strictly required to get the code, but is required to properly 
compile the code. This is due to CMacIonize's strict reproducibility 
checks that ensure all production runs use a clean version of the code.

CMacIonize supports various input and output formats, but the most 
commonly used format is [HDF5](https://www.hdfgroup.org/), a binary 
format that is optimised for massively parallel systems and that 
interfaces nicely with Python through the [h5py](https://www.h5py.org/) 
module. The code can be compiled without HDF5, but then a large number 
of input and output modules will not be available. It is therefore 
highly recommended that you install HDF5.

On Linux, you can install HDF5 using your package manager (choose the 
development package, e.g. `libhdf5-dev` or equivalent). On Mac OSX, you 
need to manually download, configure and compile the library from the 
HDF5 website. More detailed Mac OSX instructions are available in the 
[CMacIonize 
README](https://github.com/bwvdnbro/CMacIonize#dependencies).

# Obtaining the code

You can obtain CMacIonize in various ways, depending on how you want to 
use the code. The easiest way is to download a `.zip` or `.tar.gz` 
archive containing the latest stable version of the code from 
<https://github.com/bwvdnbro/CMacIonize/releases>. Unzip or untar this 
archive in a location of choice. We will call `location/CMacIonize` the 
*root folder* in what follows.

If you want a more recent version of the code, you can also `clone` the 
repository using the `git` command. Execute the following command from your
location of choice for an *anonymous* download:

```
git clone https://github.com/bwvdnbro/CMacIonize.git
```

If you plan to contribute to the code, it is advisable to create an 
account on [github](https://github.com), and then *fork* the repository 
on <https://github.com/bwvdnbro/CMacIonize>. You can then clone the 
repository using

```
git clone git@github.com:<account name>/CMacIonize.git
```

Where `<account name>` is your github account name. This will allow you 
to keep track of your changes in your version of the online repository, 
and will allow you to submit *pull requests* to include your changes 
into the parent repository.

Once you have *cloned* the repository, you should have a folder 
`location/CMacIonize`, which we will again call the *root folder*.

By default, git clones the `master` branch of the repository. This 
branch is generally the same as the latest stable release of the code. 
To get a more recent version of the code that contains unreleased new 
features, you can switch to the `stable` branch. Move to the root folder 
and execute the following command:

```
git checkout stable
```

The `stable` branch is generally ahead of the `master` branch, but 
should still be stable, i.e. not contain any untested features. An even 
more recent code version is contained in the `development` branch. This 
is the branch used by the main code developers, and is not guaranteed to 
always work. Other branches are used by individual developers during 
development and should not be used.

# Configuring the code

Once you have downloaded or cloned the code, you should have the code 
files present in the *root folder* (see above). Now you need to 
configure the code, i.e. tell the code which compiler and which version 
of the required libraries you want to use, and which specific additional 
features you want to activate. This is done using CMake. At the end of 
the configuration process, CMake will generate a `Makefile` that can 
then be used to compile the code.

To run CMake, you first need to set up a *build directory*. This could 
in principle be any folder in your system (including the root folder), 
but we advise you to use a sub-folder of the root folder, e.g. 
`root/build`. Once you have set up a build folder, open it on the 
command line, and run the following command:

```
cmake <path to root folder>
```

CMake will automatically detect the `CMakeLists.txt` configuration file 
and use it to set up the CMacIonize build environment. The output of the 
CMake command will tell you if something went wrong during configuration 
and will suggest possible solutions.

By default, CMake will use the default system compilers and libraries, but you
can instruct it to use different versions. To e.g. configure your code to use
the Intel compiler rather than the default GCC compiler on a Linux system, you
can run

```
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc <path to root folder>
```

Note that you can generally not change the compiler once you have 
already run CMake, because CMake saves it in a local cache.

Apart from changing compiler and library versions, CMake also allows you 
to activate specific features in the code, or to change code 
optimization levels. All these features will be covered in a future 
advanced tutorial.

# Compilation

Once you successfully configured the code, the *build folder* will 
contain a `Makefile`. To compile the code, you can run

```
make
```

To compile with multiple threads in parallel, run

```
make -j <number of parallel threads>
```

This speeds up the compilation process quite significantly.

If compilation succeeded (`make` will tell you if something went wrong), 
the `CMacIonize` executable will be created in `build folder/rundir/`.

CMacIonize ships with a large set of unit tests that check the accuracy of
the code. To compile and run these, run

```
make check
```

Additional `make` options will be detailed in a future advanced 
tuturial.

# Example run

At the end of the compilation process, you should have the executable 
`build folder/rundir/CMacIonize`. To start using it, you can run one of 
the example benchmark problems, which are automatically set up in `build 
folder/rundir/benchmarks/` by CMake. As a first run, we suggest running 
the Strömgren test, which consists of a single ionising source that 
ionises a Strömgren bubble inside a uniform box. Note that this test 
only works if you compiled the code with HDF5 support (none of the 
benchmark tests will work without HDF5).

Open `build folder/rundir/benchmarks/stromgren/` from the command line. 
In a clean CMacIonize build, this folder should contain three files: the 
parameter file `stromgren.param`, an analysis script `stromgren.py`, and 
a text file `stromgren.txt` containing some more information about the 
test. Assuming the CMacIonize executable is in its default location, you 
can run the benchmark test using

```
../../CMacIonize --params stromgren.param --threads <number of parallel threads>
```

You can omit the `--threads` parameter to run using a single thread. 
This test should take under a minute using 4 parallel threads.

At the end of the run, a number of additional files should have been 
created: two snapshot files `stromgren_000.hdf5` and 
`stromgren_020.hdf5`, and a file `parameters-usedvalues.param` that 
contains the parameter values that were actually used by the code 
(including parameters for which default values were used).

To check that the run produced the expected result, run the analysis 
script `stromgren.py` (this requires Python, `numpy`, `scipy`, 
`matplotlib` and `h5py`). This script will generate a radial neutral 
fraction profile for both snapshot files and also plot the expected 
radius of the Strömgren sphere. If all went well, the second snapshot 
should reproduce the expected radius reasonably well.

# Further reading

Now that you have set up the code and run your first example, you can 
start using the code for more advanced applications. Right now, there 
are no additional tutorials to help you with this, but these will appear 
in the future.
