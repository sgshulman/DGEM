# DGEM

![Build Ubuntu](https://github.com/sgshulman/DGEM/workflows/Build%20Ubuntu/badge.svg?branch=master&event=push)

A simple three-dimensional dust continuum radiative transfer code demonstrating directions grid enumeration method advantages.

This code is based on mcpolar program by Kenneth Wood (http://www-star.st-and.ac.uk/~kw25/research/montecarlo/montecarlo.html). The code was translated from FORTRAN to C++, refactored and improved. Now it realizes two radiative transfer techniques: Monte Carlo method and directions grid enumeration method (DGEM). DGEM uses precalculated directions of the photons propagation instead of the random ones to speed up the calculations process.

The physics and method detail will be presented in a paper which is in preparation now.

The provided makefile allows compiling the program on Linux.
Cmake configuration can be used to build the program and for IDEs.

PLOT3.plt is a gnuplot script for plotting the resulting images.

A utility for results comparison is differ. It is in a directory differ. Compiled differ utility requires two filenames as command arguments. The difference is called "refsum".

## DGEM configuration
The program use parameters.json configuration file with following sections:

#### method_parameters

Common settings
- fMonteCarlo &mdash; is True or False. Set False for using DGEM or set True for using Monte Carlo method
- taumin &mdash; is the minimum optical depth which gives scatterings
- nscat &mdash; number of scatterings considered

Monte Carlo parameters
- nphotons &mdash; total number of photon packets for Monte Carlo method
- iseed &mdash; initialization parameter for random number generator in Monte Carlo method

DGEM parameters
- PrimaryDirectionsLevel &mdash; primary directions number is 5 * 4 ^ PrimaryDirectionsLevel
- SecondaryDirectionsLevel &mdash; primary directions number is 5 * 4 ^ SecondaryDirectionsLevel
- NumOfPrimaryScatterings &mdash; number of scatterings in every direction in the primary grid
- NumOfSecondaryScatterings &mdash; number of scatterings in every direction in secondary grid
- MonteCarloStart &mdash; number of scatterings of the photon package after which Monte Carlo
                            method is used instead of DGEM (1 is recommended)

#### dust
- kappa &mdash; The extinction opacity
- albedo &mdash; Single scattering albedo
- hgg &mdash; Henyey-Greenstein phase function anisotropy
- pl &mdash; Peak linear polarization
- pc &mdash; Peak circular polarization
- sc &mdash; Skew Factor
      
### stars
Is a list of sources with 4 parameters:
- x, y, z &mdash; source coordinates
- l &mdash; source luminosity

Physical and disk parameters are discussed in the paper.

## Used third-party library

* **[nlohmann/json](https://github.com/nlohmann/json)** to parse a configuration file
