---
layout: tutorial
title: Computing line emissivities
date: 2019-12-04
level: basic
---

{% assign angstrom = "&#8491;" %}
{% assign micron = "&mu;m" %}

{% include toc.html %}

When run in `--emission` mode, CMacIonize can be used to compute line 
emissivities for a snapshot from a post-processing run. This assumes 
that you have already run a full post-processing Monte Carlo 
photoionization simulation and that you have a snapshot that contains 
accurate temperatures and ionic fractions for all cells.

In this mode, CMacIonize requires two additional command line arguments:
 * `--params`" The name of the parameter file, as usual.
 * `--file`: The name of the snapshot file.

The parameter file will have a completely different contents from the 
usual parameter file, as explained below. If the calculation of the 
emissivities is successful, one or more new datasets containing the 
result will be appended to the snapshot file.

# Parameter file

Unlike the usual parameter file, the emission-mode parameter file only 
has one parameter block, called `EmissivityValues`. Currently, the 
following parameters are available:

```
EmissivityValues:
  BaHigh: false # (default value)
  BaLow: false # (default value)
  CIII_1908: false # (default value)
  CII_158mu: false # (default value)
  CII_2325: false # (default value)
  HII: false # (default value)
  Halpha: false # (default value)
  Hbeta: false # (default value)
  HeI_5876: false # (default value)
  Hrec_s: false # (default value)
  NIII_57mu: false # (default value)
  NII_122mu: false # (default value)
  NII_5755: false # (default value)
  NII_6548: false # (default value)
  NII_6584: false # (default value)
  NeIII_15mu: false # (default value)
  NeIII_3869: false # (default value)
  NeIII_3968: false # (default value)
  NeII_12mu: false # (default value)
  OIII_4363: false # (default value)
  OIII_4959: false # (default value)
  OIII_5007: false # (default value)
  OIII_52mu: false # (default value)
  OIII_88mu: false # (default value)
  OII_3727: false # (default value)
  OII_7325: false # (default value)
  OI_6300: false # (default value)
  SIII_19mu: false # (default value)
  SIII_33mu: false # (default value)
  SIII_6213: false # (default value)
  SIII_9405: false # (default value)
  SII_4072: false # (default value)
  SII_6725: false # (default value)
  SIV_10mu: false # (default value)
  WFC2_F439W: false # (default value)
  WFC2_F555W: false # (default value)
  WFC2_F675W: false # (default value)
  avg_T: false # (default value)
  avg_T_count: false # (default value)
  avg_nH_nHe: false # (default value)
  avg_nH_nHe_count: false # (default value)
```

As per usual, CMacIonize will read the parameter file and then generate 
a file called `ORIGINAL_NAME.used-values` (with `ORIGINAL_NAME` the name 
of the input parameter file; including extension) that contains the 
parameters that were actually used (including all default values).

All parameters in the file correspond to a specific line emissivity 
diagnostic, and the value of the parameter decides whether or not this 
line will be computed. By default, no lines are computed.

Below is an overview of the line diagnostics currently supported:

## `BaHigh`, `BaLow`

Legacy values from Kenny Wood's code, currently undocumented.

## `CIII_1908`

This corresponds to the 1908 {{ angstrom }} transition line from the 
(degenerate) first excited state of double ionised carbon to its ground 
state (see Osterbrock & Ferland, 2006; table 3.8, 3P0:3P1:3P2 to 1S0).

## `CII_158mu`, `CII_2325`

These correspond to the 158 {{ micron }} hyperfine transition line 
between the degenerate levels in the ground state of ionised carbon 
(2P3/2 to 2P1/2), and the 2325 {{ angstrom }} transition from the 
(degenerate) first excited state to the (degenerate) ground state (see 
Osterbrock & Ferland, 2006; table 3.9, 4P1/2:4P3/2:4P5/2 to 
2P1/2:2P3/2).

## `HII`

Legacy values from Kenny Wood's code, currently undocumented.

## `Halpha`

H&alpha; emissivity, based on the power law fit provided by Osterbrock & 
Ferland (2006), table 4.1.

## `Hbeta`

H&beta; emissivity, based on a power law fit to the Storey & Hummer 
(1995) data.

## `HeI_5876`

5876 {{ angstrom }} emission line from neutral helium, based on a power 
law fit provided by Osterbock & Ferland (2006), table 4.6.

## `Hrec_s`

Hydrogen recombination rate, as a number of recombining ions per second, 
based on the fitting formula in Verner & Ferland (1996), table 1.

## `NIII_57mu`

Corresponds to the 57 {{ micron }} transition from the excited state to 
the ground state in double ionised nitrogen (Osterbrock & Ferland, 2006; 
table 3.9, 2P3/2 to 2P1/2).

## `NII_122mu`, `NII_5755`, `NII_6548`, `NII_6584`

These correspond to transitions between excited states of ionised nitrogen
(Osterbrock & Ferland, 2006; table 3.12):
 * the 122 {{ micron }} hyperfine transition between the degenerate 
levels of the first excited state (3P2 to 3P1)
 * the 5755 {{ angstrom }} transition between the third and second 
excited state (1S0 to 1D2)
 * the 6548 {{ angstrom }} transition between the second excited state 
and the lower level of the (degenerate) first excited state (1D2 to 3P1)
 * the 6584 {{ angstrom }} transition between the second excited state 
and the higher level of the (degenerate) first excited state (1D2 to 
3P2)

## `NeIII_15mu`, `NeIII_3869`, `NeIII_3968`

These correspond to transitions between excited states of double ionised 
neon (Osterbrock & Ferland, 2006; table 3.14):
 * the 15 {{ micron }} hyperfine transition between the first and second 
level in the degenerate ground state (3P1 to 3P2)
 * the 3869 {{ angstrom }} transition between first excited state and 
the lowest level of the (degenerate) ground state (1D2 to 3P2)
 * the 3968 {{ angstrom }} transition between the first excited state 
and the second level of the (denegerate) ground state (1D2 to 3P1)

## `NeII_12mu`

12 {{ micron }} hyperfine transition between the two levels in the 
ground state of ionised neon (Osterbrock & Ferland, 2006; table 3.11, 
2P1/2 to 2P3/2).

## `OIII_4363`, `OIII_4959`, `OIII_5007`, `OIII_52mu`, `OIII_88mu`

These correspond to transitions between excited states of double ionised 
oxygen (Osterbrock & Ferland, 2006; table 3.12):
 * the 4363 {{ angstrom }} transition between the second and first 
excited state (1S0 to 1D2)
 * the 4959 {{ angstrom }} transition between the first excited state 
and the second level of the (degenerate) ground state (1D2 to 3P1)
 * the 5007 {{ angstrom }} transition between the first excited state 
and the highest level of the (degenerate) ground state (1D2 to 3P2)
 * the 52 {{ micron }} hyperfine transition between the third and second 
level in the degenerate ground state (3P2 to 3P1)
 * the 88 {{ micron }} hyperfine transition between the second and first 
level in the degenerate ground state (3P1 to 3P0)

## `OII_3727`, `OII_7325`

These correspond to transitions between excited state of ionised oxygen 
(Osterbrock & Ferland, 2006; table 3.13):
 * the sum of the 3726 and 3728.8 {{ angstrom }} transition lines 
between the (degenerate) first excited state and the ground level 
(2D5/2:2D3/2 to 4S3/2)
 * the sum of the 7319.9, 7330.7, 7318.8 and 7329.6 {{ angstrom }} 
transition lines between the (degenerate) second and first excited 
states (2P3/2:2P1/2 to 2D5/2:2D3/2)

## `OI_6300`

Sum of the 6300.3 and 6363.8 {{ angstrom }} transition lines between the 
first excited level of neutral oxygen and its (degenerate) ground state 
(Osterbrock & Ferland, 2006; table 3.14, 1D2 to 3P2:3P1).

## `SIII_19mu`, `SIII_33mu`, `SIII_6213`, `SIII_9405`

These correspond to transitions between excited levels of double ionised 
sulphur (Osterbrock & Ferland, 2006; table 3.12):
 * the 19 {{ micron }} hyperfine transition between the third and second 
level of the (degenerate) ground state (3P2 to 3P1)
 * the 33 {{ micron }} hyperfine transition between the second and first 
level of the (degenerate) ground state (3P1 to 3P0)
 * the 6213 {{ angstrom }} transition between second and first excited 
state (1S0 to 1D2)
 * the sum of the 9531 and 9068.9 {{ angstrom }} transitions between the 
first excited state and the (degenerate) ground state (1D2 to 2P1:3P2)

## `SII_4072`, `SII_6725`

These correspond to transitions between excited levels of ionised 
sulphur (Osterbrock & Ferland, 2006; table 3.13):
 * the sum of the 4068.6 and 4076.4 {{ angstrom }} transitions between 
the (degenerate) second excited level and the ground state (2P1/1:2P3/2 
to 4S3/2)
 * the 6725 {{ angstrom }} transition between the (degenerate) first 
excited level and the ground state (2D3/2:2D5/2 to 4S3/2)

## `SIV_10mu`

Transition between the first excited state and the ground state in 
triple ionised sulphur (Osterbrock & Ferland, 2016; table 3.10; 2P3/2 to 
2P1/2).

## `WFC2_F439W`, `WFC2_F555W`, `WFC2_F675W`

Sums of all the lines that fall within the wavelength ranges of the HST 
WFC2 filters F439W, F555W and F675W, useful to make synthetic HST 
images.

## `avg_T`, `avg_T_count`

Variables required to compute the density weighted average temperature 
of the ionised particles in each cell.

## `avg_nH_nHe`, `avg_nH_nHe_count`

Variables required to compute the average product of the hydrogen and 
helium ionised densities.
