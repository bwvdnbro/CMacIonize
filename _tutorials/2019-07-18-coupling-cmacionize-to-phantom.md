---
layout: tutorial
title: Coupling CMacIonize to PHANTOM
date: 2019-07-18
level: advanced
---

{% include toc.html %}

[PHANTOM](https://phantomsph.bitbucket.io/) is a public SPH code that 
can be coupled to CMacIonize to perform RHD simulations. This setup uses 
the CMacIonize Fortran library. At specific times during the SPH 
calculation, the particle positions, masses and smoothing lengths are 
passed on to the library. CMacIonize then constructs a Voronoi grid 
using the particle positions as Voronoi cell generators and maps the 
particle densities onto this grid. This grid is then used for a 
photoionization simulation. At the end of the simulation, the resulting 
hydrogen neutral fractions are returned to PHANTOM and used to update 
the thermal properties of the SPH particles.

While this sounds all very good, there are a few technical challenges 
when setting up these simulations, ranging from installing both PHANTOM 
and CMacIonize, over making sure both codes can talk to each other, to 
choosing good values for the many parameters that control the hybrid 
code. This tutorial aims to address these issues.

# Installation

## Obtaining and installing PHANTOM

A very comprehensive manual on getting started with PHANTOM can be found 
on <https://bitbucket.org/danielprice/phantom/wiki/getstarted>. Below 
are some known pitfalls:
 - When checking out the latest version of the code from the BitBucket 
repository, make sure to check that you have a stable version (this can 
be verified from the git tag using `git describe`). The safest way to do 
this is by using the `stable` branch of the repository (`git checkout 
stable`).
 - Make sure to set the correct environment variables before invoking 
any `make` command or running `phantom`: the `OMP_SCHEDULE` and 
`OMP_STACKSIZE` variables described at the start of the manual and the 
`SYSTEM` variable that specifies the Fortran compiler to use.
 - It is good practice to keep the PHANTOM source directory clean; put 
your simulation files in another directory (e.g. `phantom-build`).

## Obtaining and installing CMacIonize

This is already covered in [another tutorial]({% link 
_tutorials/2019-01-11-getting-started-with-cmacionize.md %}). Note that 
the main development of the coupling module is done in a 
[fork](https://github.com/mapetkova/CMacIonize) of the main repository. 
So it might be worthwhile to use that version of the code. Since a fork 
acts as a completely independent repository, this means you will need to 
clone it into a separate folder outside the main CMacIonize folder.

## Telling PHANTOM where to find CMacIonize

More information here.

# Running a PHANTOM+CMacIonize RHD simulation

More information here.
