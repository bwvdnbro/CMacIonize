---
layout: tutorial
title: Creating your own parameter files
date: 2019-01-14
level: basic
---

{% include toc.html %}

The *parameter file* is the most important component of any CMacIonize 
simulation, as it contains all the parameters that determine what your 
simulation will do and how it will do this, including
 * input parameters: what is the initial condition of the simulation?
 * run time parameters: how do you want to run your simulation?
 * output parameters: what do you want to output and how?

The parameter file is written in `YAML` format (<https://yaml.org/>), 
and consists of a (large) number of *blocks*. Each block maps to a 
component of the simulation (and usually also a specific *module* in the 
code). The majority of the blocks have a *field* called `type`, which 
allows you to choose between different variations (*implementations*) of 
the same module, e.g. different initial conditions or different ISM 
properties. Furthermore, the parameter file supports units (in fact, you 
are *forced* to use units for all physical quantities that have a unit).

In this tutorial, you will learn how to set up your own parameter file. 
*It will not cover all available types for each module, as this list 
constantly grows.* Instead, it will show you how to find an up-to-date 
list of all types, and how you can use the output of CMacIonize to 
fine-tune your parameter file.

# Default parameter file

In order to run CMacIonize, you need to pass a parameter file to the 
program when running it, i.e.

```
./CMacIonize --params <name of parameter file>
```

This does however not mean that the parameter file you provide has to 
contain any content. CMacIonize has a set of default parameters for each 
variable in the parameter file, and it will use these if it cannot find 
that specific variable (that is what *default* means). You can use this 
fact to your advantage to make CMacIonize output a default parameter 
file, that contains all available variables and their default values.

First, create an (empty) parameter file. On Linux systems, this can be 
done using `touch`:

```
touch empty.param
```

Now run CMacIonize using this parameter file, and force it to stop after 
the program initialisation by using the `dry-run` option:

```
./CMacIonize --params empty.param --dry-run
```

This will start a simulation, read in all necessary parameters and 
allocate all variables in memory, and then stop the program. This is 
perfect to test whether your simulation actually fits on a specific 
machine, and to check your parameters. After CMacIonize has read all 
parameters, it will write a file called `parameters-usedvalues.param`, 
that contains all parameters that were read during initialisation, and 
the values that were assigned to these variables.

At the moment of writing, this is what the default 
`parameters-usedvalues.param` file looks like:

```
# file written on 28/02/2020, 17:38:17.
AbundanceModel:
  C: 0 # (default value)
  He: 0 # (default value)
  N: 0 # (default value)
  Ne: 0 # (default value)
  O: 0 # (default value)
  S: 0 # (default value)
  type: FixedValue # (default value)
ContinuousPhotonSource:
  type: None # (default value)
ContinuousPhotonSourceSpectrum:
  frequency: 3.28847e+15 Hz # (default value)
  total flux: -1 m^-2 s^-1 # (default value)
  type: Monochromatic # (default value)
CrossSections:
  type: Verner # (default value)
DensityFunction:
  density: 1e+08 m^-3 # (default value)
  neutral fraction H: 1e-06 # (default value)
  temperature: 8000 K # (default value)
  type: Homogeneous # (default value)
DensityGrid:
  number of cells: [64, 64, 64] # (default value)
  type: Cartesian # (default value)
DensityGridWriter:
  padding: 3 # (default value)
  prefix: snapshot # (default value)
  type: Gadget # (default value)
DensityGridWriterFields:
  Coordinates: 1 # (default value)
  CosmicRayFactor: 0 # (default value)
  NeutralFractionC+: 0 # (default value)
  NeutralFractionC++: 0 # (default value)
  NeutralFractionH: 1 # (default value)
  NeutralFractionHe: 0 # (default value)
  NeutralFractionN: 0 # (default value)
  NeutralFractionN+: 0 # (default value)
  NeutralFractionN++: 0 # (default value)
  NeutralFractionNe: 0 # (default value)
  NeutralFractionNe+: 0 # (default value)
  NeutralFractionO: 0 # (default value)
  NeutralFractionO+: 0 # (default value)
  NeutralFractionS+: 0 # (default value)
  NeutralFractionS++: 0 # (default value)
  NeutralFractionS+++: 0 # (default value)
  NumberDensity: 1 # (default value)
  Temperature: 0 # (default value)
DensityMask:
  type: None # (default value)
DiffuseReemissionHandler:
  type: None # (default value)
IonizationSimulation:
  enable trackers: false # (default value)
  number of iterations: 10 # (default value)
  number of photons: 100000 # (default value)
  number of photons first loop: 100000 # (default value)
  output folder: . # (default value)
  random seed: 42 # (default value)
PhotonSourceDistribution:
  luminosity: 4.26e+49 Hz # (default value)
  position: [0 m, 0 m, 0 m] # (default value)
  type: SingleStar # (default value)
PhotonSourceSpectrum:
  frequency: 3.28847e+15 Hz # (default value)
  total flux: -1 m^-2 s^-1 # (default value)
  type: Monochromatic # (default value)
RecombinationRates:
  type: Verner # (default value)
SimulationBox:
  anchor: [-1.543e+17 m, -1.543e+17 m, -1.543e+17 m] # (default value)
  periodicity: [false, false, false] # (default value)
  sides: [3.086e+17 m, 3.086e+17 m, 3.086e+17 m] # (default value)
TemperatureCalculator:
  PAH heating factor: 0 # (default value)
  cosmic ray heating factor: 0 # (default value)
  cosmic ray heating limit: 0.75 # (default value)
  cosmic ray heating scale length: 4.11466e+19 m # (default value)
  do temperature calculation: false # (default value)
  epsilon convergence: 0.001 # (default value)
  maximum number of iterations: 100 # (default value)
  minimum ionized temperature: 4000 K # (default value)
  minimum number of iterations: 3 # (default value)
```

These values correspond to a variation on the Strömgren benchmark test 
that is part of the benchmarks. The blocks and values are ordered 
alphabetically and are hence not necessarily in a sensible order. Each 
variable has two values: the actually used value, which for physical 
quantities is always given in SI units, and the value that was provided 
in the parameter file (or *default value* if no value was given in the 
parameter file; this is true for all variables in this example). If a 
value was not used this will also be clearly indicated. When a parameter 
specifies the name of a file in the local file system (e.g. when using a 
module that reads a file to set up the initial conditions of the 
simulation), the `parameters-usedvalues.param` will also include a 
checksum of that file. If the checksum is provided, the code will throw 
an error if that checksum does not match the checksum for the file with 
that name that is present.

You can now use this default parameter file to start creating your own 
parameter file: copy the file (e.g. `cp parameters-usedvalues.param 
my_params.param`) to make sure it is not overwritten the next time you 
run CMacIonize, and edit it to your liking. Then re-run CMacIonize and 
inspect the new `parameters-usedvalues.param` to check if your new 
parameters were read correctly (or to get a new list of default values 
if you changed a block *type*). Keep doing this until your are satisfied 
with the result.

# Parameter blocks

Below is a full description of the various blocks in the default 
parameter file, to give you a better idea of what the different blocks 
and parameters do. Note that new parameter blocks are added to the code 
regularly, so that this overview is not necessarily up-to-date. New 
blocks are always created with default values that guarantee that the 
Strömgren benchmark is the default simulation, so if you are unsure 
about a parameter block, it is usually fine if you just keep the default 
values.

## AbundanceModel

This block specifies the model to use for the abundances of all elements 
that CMacIonize can keep track of, as the corresponding number density 
normalised to the number density of hydrogen. These are used by various 
modules throughout the code and affect the cooling rates used for 
temperature calculations, as well as the optical depth calculation 
during photon packet propagation. The default type, `FixedValue`, allows 
to set the abundances of all elements manually. For RHD simulations, it 
is a good idea to set all these values to 0, resulting in a 
hydrogen-only ISM.

## ContinuousPhotonSource/ContinousPhotonSourceSpectrum

These blocks allow control over non-discrete photon sources, i.e. smooth 
luminosity distributions or external radiation fields. By default, these 
sources are disabled. Available `ContinuousPhotonSourceSpectrum` types 
are identical to those for `PhotonSourceSpectrum`.

## CrossSections

This block sets the photoionisation cross sections for all ions tracked 
by the code. The default `Verner` type corresponds to a self-consistent 
set of values from literature. The alternative `FixedValue` type allows 
you to manually set the cross section for the individual ions.

## DensityFunction

This parameter controls the way grid cell values are initialised at the 
start of the simulation and hence sets the initial condition for the 
simulation. All types work in the same way: they create a *function* 
that takes a cell of the grid as input and that returns the initial 
values for the number density, temperature, velocity and neutral 
fractions of that cell. The default type returns constant values for 
each of these quantities (and a velocity that is always 0). Alternatives 
include specific density profile models (e.g. 
`DiscPatchDensityFunction`) or density profiles read from other files 
(e.g. `FLASHSnapshotDensityFunction`).

This is one of the blocks you are most likely to change.

## DensityGrid

This block controls the grid. The type selects the type of grid (at the 
time of writing, the available options are `Cartesian` for a regular 
grid, `AMR` for an adaptive grid, and `Voronoi` for an unstructured 
grid). All types have parameters that control the number of grid cells, 
for the default Cartesian grid this is an array that contains the number 
of cells in the 3 coordinate directions (the total number of cells is 
the product of this number).

## DensityGridWriter

This block controls the output of the simulation. The type parameter 
selects one of the available output modules (`Gadget` for the default 
Gadget2-style HDF5-output, `AsciiFile` for a simple text file). The 
`prefix` and `padding` parameters control respectively the name of the 
snapshot files and the number of digits in the automatic counter for 
these files.

## DensityGridWriterFields

This block also concerns output and controls the variables that are 
actually written to every output file. Note that at the moment of 
writing, this feature is only fully supported by the default 
`GadgetDensityGridWriter`. The values for the various parameters are 
booleans; `true/yes/on/y/1` (and variants with upper case letters) will 
result in a block with the corresponding name being written to the 
output file, `false/no/off/n/0` has the reverse effect.

Note that the available fields depend on the type of simulation. RHD 
simulations have additional fields for hydrodynamical variables that are 
not present for post-processing simulations.

## DensityMask

This block makes it possible to apply an additional filter to the 
density field generated by the `DensityFunction`, e.g. redistribute the 
density values to generate a clumpy, fractal substructure 
(`FractalDensityMask`). By default, this is disabled. Note that this 
feature might be deprecated in the future.

## DiffuseReemissionHandler

This parameter block determines how diffuse reemission by the ISM is 
modeled, i.e. if and how ionizing radiation that was absorbed can be 
reemitted as ionizing radiation. The default type, `None`, disables 
diffuse reemission. Type `Physical` is a good choice for simulations 
that do need to include diffuse reemission.

## IonizationSimulation

This block contains the main run time parameters for the simulation, and 
only exists for a post-processing simulation. For an RHD simulation, it 
is replaced by two other blocks, called 
`RadiationHydrodynamicsSimulation` and `HydroIntegrator`.

The parameters listed in this block (and equivalently in 
`RadiationHydrodynamicsSimulation`) control the number of iterations and 
the number of photon packets for the Monte Carlo algorithm, the random 
seed and the location where all output is stored on the system. For 
`RadiationHydrodynamicsSimulation`, there are additional parameters that 
control the time line of the simulation and the frequency at which the 
radiation transfer is performed or output is written.

The `HydroIntegrator` block (absent in this example) controls the 
parameters of the hydrodynamics solver for an RHD simulation.

## PhotonSourceDistribution

This block controls the distribution of discrete sources of radiation 
and can hence also be considered part of the initial condition. As for 
the `DensityFunction`, all types work as a sort of *function* that 
generates a number of discrete source positions, and that also outputs 
the total luminosity and weight of each source. The code treats all 
these sources in the same way and is set up so that it can deal with 
them in a general way. The default type is a single source of radiation. 
Other available types generate source positions randomly from an 
analytic distribution function (e.g. `SILCCPhotonSourceDistribution`) or 
read them from an input file (e.g. 
`GadgetSnapshotPhotonSourceDistribution`).

Again, you probably need to adapt this block to the specific simulation 
you want to run.

## PhotonSourceSpectrum

This block controls the spectrum of the ionising sources. Among the 
alternatives for the default monochromatic spectrum are a black-body 
spectrum (`PlanckPhotonSourceSpectrum`) and a model stellar spectrum 
(`WMBasicPhotonSourceSpectrum`).

## RecombinationRates

This block sets the recombination rates for the various ions tracked by 
the code. The default `VernerRecombinationRates` contains 
self-consistent recombination rate data from literature; the alternative 
`FixedValueRecombinationRates` allows you to manually specify the 
recombination rate for individual ions.

## SimulationBox

This block sets the dimesnions of the simulation box. The `anchor` 
parameter specifies the bottom front left corner of the box (lowest `x`, 
`y` and `z` coordinates), the `sides` parameter specifies the extent of 
the box in the three coordinate directions. The `periodicity` parameter 
specifies what happens to photon packets that leave the box in each of 
the three coordinate directions: a value of `true` (or `yes/y/on/1`, see 
above) means these packets reenter the box on the other side, as if the 
box was surrounded by an infinite amount of identical copies in that 
coordinate direction. Note that this has a significant impact on the run 
time of the simulation.

Also note that the `periodicity` parameter in this block only affects 
the Monte Carlo algorithm. Boundary conditions for the hydrodynamics 
solver in RHD simulations are set in the `HydroIntegrator` block. The 
reason for this is that there are a lot of additional options for 
hydrodynamical boundaries.

## TemperatureCalculator

This last block is used for the full mode of the Monte Carlo 
photoionisation algorithm in which temperatures are updated 
self-consistently by keeping track of photoheating and recombination 
line cooling. The `do temperature calculation` parameter controls 
whether or not this full mode is active (it is disabled by default). 
Additional parameters control the details of the recursive algorithm and 
additional heating sources.

Note that enabling the full mode of the code has a significant impact on 
the run time of the code. Furthermore, a larger number of photon packets 
is required to accurately model the various coolants than is required to 
just get an accurate ionisation structure for hydrogen and helium.

# Additional documentation

All blocks are extensively documented in the [doxygen 
documentation](https://users.ugent.be/~bwvdnbro/CMacIonize/) of the C++ 
class with the same name. To get a complete and up-to-date list of all 
available types for a specific block, look at the documentation for the 
corresponding *factory* class (e.g. `DensityFunctionFactory` for a list 
of all available `DensityFunction` types).

Note that the online doxygen documentation is not guaranteed to be 
up-to-date. To get up-to-date documentation for the code version you are 
using, manually compile the documentation using `make doc` (this 
requires a working `doxygen` installation).
