---
layout: tutorial
title: Command line parameters
date: 2019-01-15
level: basic
---

{% include toc.html %}

CMacIonize has a number of command line parameters that control the 
execution of the program. These command line parameters are not 
equivalent to the parameters in the parameter file, as changing them 
does not alter the result of the simulation, but only affects the way 
CMacIonize runs. They allow you to change the resources that CMacIonize 
uses, control the on-screen output, or select specific modes of the 
code. This tutorial gives an overview that is complete at the time of 
writing.

# Compulsory command line parameters

CMacIonize has only one compulsory command line parameter:

```
--params <name of parameter file>
```

With this parameter, you can specify the name of the parameter file. 
This can be given as an absolute path, or as a relative path (in which 
case the directory that acts as working directory for the `CMacIonize` 
program call will be used as reference path).

Note that the parameter file does not need to contain any actual 
parameters ([see the parameter file tutorial]({% link 
_tutorials/2019-01-14-creating-your-own-parameter-files.md %})).

# (Optional) program mode command line parameters

By default, the code will run in post-processing mode and perform a 
single photoionisation loop to compute the ionisation structure and 
optionally the temperature of the given ISM density structure. Two other 
modes are available, they are activated by specifying appropriate 
command line parameters:

```
--rhd
```

This parameter selects the radiation hydrodynamics (RHD) mode of the 
code. The code will perform a time integration of the initial ISM 
density structure using a combination of hydrodynamical integration 
steps and intermediate photoionisation steps.

Note that this mode (de-)activates a number of parameter blocks in the 
parameter file, so generally a parameter file that works for a 
post-processing simulation will not work for an RHD simulation.

```
--dusty-radiative-transfer
```

This mode of the code has a very specific purpose, as it activates an 
algorithm that performs a post-processing simulation of dust scattering 
through the given ISM structure. This mode is only present because this 
specific use case was required for teaching purposes, and should be 
considered as an unstable hack. This mode only accepts a subset of the 
available modules and is not suited for scientific simulations.

```
--emission
```

_NEW_ This mode of the code computes line emission strengths for a 
snapshot that was generated using the full version of the code in 
post-processing, and that contains accurate temperature values and ionic 
fractions for all elements.

When run in this mode, CMacIonize will simply read the snapshot that is 
passed on to the program and compute the requested line emissivities 
using the temperature and ionic fractions present in the snapshot file. 
It will then append the emissivities as additional datasets to the 
snapshot.

## Command line parameters that affect performance

```
--threads <number of threads to use>
```

This parameter sets the number of shared-memory parallel threads that 
will be used by the code. Both the Monte Carlo algorithm and the 
hydrodynamics solver scale well on a shared-memory system, so setting 
this value to a sufficiently high number will significantly decrease the 
run time of your simulation.

Note that the number of threads is capped at the number of available 
cores on the system, so selecting a higher value will not work. This 
reference value is determined during the code configuration, so always 
make sure to configure and compile the code on the same system as you 
will be using for production runs!

```
--dry-run
```

This parameter can be used to perform a *dry run* of the simulation. The 
code will read all the parameters and set up the entire model, and will 
then gracefully exit without actually running the simulation. This is 
ideal to test whether a simulation can actually run on a specific system 
(e.g. to check if there is enough memory), to check if all parameters 
are set to sensible and consistent values, and to iteratively set up a 
parameter file as explained in [the parameter file tutorial]({% link 
_tutorials/2019-01-14-creating-your-own-parameter-files.md %}).

## Command line parameters that affect reproducibility

```
--use-version <tag of a specific git version>
```

This parameter allows you to force the use of a specific version of the 
code. The code will check the tag returned by the `git describe --tags 
--dirty` command during its compilation against the given value, and 
will exit if the two values do not match. Note that the code will not 
attempt to actually set up the requested version. This flag is useful as 
a sanity check if you want to make sure that you are not accidentally 
using the wrong code version.

```
--dirty
```

By default, the code will only run if the tag returned by `git describe 
--tags --dirty` during its compilation does not contain the keyword 
*dirty*. In git terms, this means that the local repository does not 
contain any uncommitted changes and hence corresponds to a reproducible, 
tracked code version. We strongly advise you to stick to this default 
behaviour.

However, sometimes it is useful to force running the code, even when 
uncommitted changes are present. This is especially useful when testing 
newly developed modules or system specific bug-fixes before they are 
committed. By specifying the `--dirty` flag, you can force the code to 
still run, even when the code is *dirty*.

# Command line parameters that affect output

```
--verbose
```

By default, the code writes a fairly reasonable amount of log 
information to the screen while it is running. By using this parameter, 
you can increase the amount of info that is written.

Internally, this flag sets the *log level* (the level of logging 
information that is actually output) to the lowest possible value, 
resulting in all log messages to be written.

```
--logfile <name of file>
```

By default, log output is written to the system `stdout` (this is 
usually the terminal window in which you run the code). By specifying 
this command line option, you can reroute the log output to the file 
with the given path and name. Note that the path can be absolute or 
relative. In the latter case, the reference path is the folder that acts 
as working directory for the `CMacIonize` program call.

```
--every-iteration-output
```

This parameter mainly affects the post-processing mode of the code. By 
default, this mode only produces two snapshot output: a snapshot of the 
initial state of the system, and a snapshot at the end of the 
photoionisation loop. By specifying this flag, you can force the code to 
write a snapshot at the end of every iteration in the photoionisation 
loop.

Note that this flag does not work well for RHD simulations, and should 
not be used in this case!

```
--output-statistics
```

This parameter forces the code to compute and output statistical 
information about the photon packets at the end of each photoionisation 
iteration. The flag only has effect when running in post-processing mode 
and does not do anything in RHD mode.

# RHD only command line parameters

The RHD mode of the code has two additional command line parameters that 
do not exit in post-processing mode:

```
--restart <name of restart folder>
```

This parameter instructs the code to restart a simulation that was 
previously dumped to a restart file. Its value should match the path to 
the folder that contains the restart file.

Note that restarting ignores the parameter file (although you still need 
to provide the `--params` argument), and will continue the previous run 
as if it never stopped. The code currently does not do any version 
checks (these will be implemented at some stage), but you can expect 
that restarting will go wrong if you use different versions. Restarting 
is generally only possible on the same machine as the original run, 
although migrating to another machine works occasionally. There is no 
need to use the same value for the `--threads` parameter when 
restarting.

Restarting is useful when running simulations on a system with a time 
limit, and can also be used to recover a simulation that crashed because 
of external reasons (hardware failure, system reboot...).

```
--output-time-unit <valid time unit>
```

This parameter allows you to specify the time unit used for time step 
information in the log output. By default, all time step information is 
written in SI units (s). The given time unit needs to be supported by 
CMacIonize. Currently supported time units are *s*, *yr*, *Myr*, *Gyr*.

# Help!

The last command line parameter is probably the most useful one:

```
--help
```

In the true POSIX tradition, this parameter instructs the code to print 
a list of all available command line parameters and a short description 
of what they do. This will give you a summary of this tutorial that will 
always contain all available command line parameters.
