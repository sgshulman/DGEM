# DGEM
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
- fMonteCarlo  is 0 or 1. Set 0 for using DGEM or set 1 for using Monte Carlo method
- taumin       is the minimum optical depth which gives scatterings
- nscat        number of scatterings considered

Monte Carlo parameters
- nphotons     total number of photon packets for Monte Carlo method
- iseed        initialization parameter for random number generator in Monte Carlo method

DGEM parameters
- PrimaryDirectionsLevel      primary directions number is 5 * 4 ^ PrimaryDirectionsLevel
- SecondaryDirectionsLevel primary directions number is 5 * 4 ^ SecondaryDirectionsLevel
- NumOfPrimaryScatterings     number of scatterings in every direction in the primary grid
- NumOfSecondaryScatterings    number of scatterings in every direction in secondary grid
- MonteCarloStart               number of scatterings of the photon package after which Monte Carlo
                            method is used instead of DGEM (1 is recommended)

#### dust
- kappa --- The extinction opacity
- albedo --- Single scattering albedo
- hgg ---- Henyey-Greenstein phase function anisotropy
- pl --- Peak linear polarization
- pc --- Peak circular polarization
- sc --- Skew Factor
      
### stars
Is a list of sources with 4 parameters:
- x, y, z --- coordinates
- l --- source luminosity

Physical and disk parameters are discussed in the paper.

## Used third-party library

* **[nlohmann/json](https://github.com/nlohmann/json)** to parse a configuration file
