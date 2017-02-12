# DGEM
A simple three-dimensional dust continuum radiative transfer code demonstrating directions grid enumeration method advantages.

This code is based on mcpolar program by Kenneth Wood (http://www-star.st-and.ac.uk/~kw25/research/montecarlo/montecarlo.html). The code was translated from FORTRAN to C++, refactored and improved. Now it realizes two radiative transfer techniques: Monte Carlo method and directions grid enumeration method (DGEM). DGEM uses precalculated directions of the photons propagation instead of the random ones to speed up the calculations process.

The physics and method detail will be presented in a paper which is in preparation now.
