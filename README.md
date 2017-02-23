# DGEM
A simple three-dimensional dust continuum radiative transfer code demonstrating directions grid enumeration method advantages.

This code is based on mcpolar program by Kenneth Wood (http://www-star.st-and.ac.uk/~kw25/research/montecarlo/montecarlo.html). The code was translated from FORTRAN to C++, refactored and improved. Now it realizes two radiative transfer techniques: Monte Carlo method and directions grid enumeration method (DGEM). DGEM uses precalculated directions of the photons propagation instead of the random ones to speed up the calculations process.

The physics and method detail will be presented in a paper which is in preparation now.

The provided makefile allows to compile the program on Linux.

PLOT3.plt is a gnuplot script for plotting the resulting images.

An utility for results comparison is differ. It is in a directory differ. Compiled differ utility requires two filename as command arguments. The difference is called "refsum".

The program use param.par file with options:
Method parameters
- fMonteCarlo  is 0 or 1. Set 0 for using DGEM or set 1 for using Monte Carlo method
- taumin       is minimum optical depth which gives scatterings
- nscat        number of scatterings considered
  
Monte Carlo Parameters
- nphotons     total number of photon packets
- iseed        initialization parameter for random number generator

DGEM Parameters
- PrimaryDirectionsLevel 		primary directions number is 5 * 4 ^ PrimaryDirectionsLevel
- SecondaryDirectionsLevel	primary directions number is 5 * 4 ^ SecondaryDirectionsLevel
- NumOfPrimaryScatterings		number of scatterings in every direction in primary grid
- NumOfSecondaryScatterings	number of scatterings in every direction in secondary grid
- MonteCarloStart				    number of scatterings of the photon package after which Monte Carlo
                            method is used instead of DGEM (1 is recommended)

Stars
  nstars     is a number of stars. After it goes nstar lines with star coordinates and luminosities. The example is given.
  star=0.0 0.0 0.0 1.0
Physical and disk parameters are discussed in the paper.
