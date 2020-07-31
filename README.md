# DGEM

![Unit Tests Ubuntu](https://github.com/sgshulman/DGEM/workflows/Unit%20Tests%20Ubuntu/badge.svg?branch=master&event=push)
![Unit Tests macOS](https://github.com/sgshulman/DGEM/workflows/Unit%20Tests%20macOS/badge.svg?branch=master&event=push)
![Unit Tests Windows](https://github.com/sgshulman/DGEM/workflows/Unit%20Tests%20Windows/badge.svg?branch=master&event=push)

A simple three-dimensional dust continuum radiative transfer code demonstrating directions grid enumeration method advantages.

- [DGEM Configuration](#dgem-configuration)
  - [Method Parameters](#method-parameters)
  - [Dust](#dust)
  - [Geometry](#geometry)
  - [Stars](#stars)
  - [Observers](#observers)
- [Used third-party libraries](#used-third-party-libraries)

This code is based on mcpolar program by Kenneth Wood (http://www-star.st-and.ac.uk/~kw25/research/montecarlo/montecarlo.html). The code was translated from FORTRAN to C++, refactored and improved. Now it realizes two radiative transfer techniques: Monte Carlo method and directions grid enumeration method (DGEM). DGEM uses precalculated directions of the photons propagation instead of the random ones to speed up the calculations process.

The physics and method detail will be presented in a paper which is in preparation now.

The provided makefile allows compiling the program on Linux.
Cmake configuration can be used to build the program and for IDEs.

PLOT3.plt is a gnuplot script for plotting the resulting images.

A utility for results comparison is differ. It is in a directory differ. Compiled differ utility requires two filenames as command arguments. The difference is called "refsum".

## DGEM configuration
The program use parameters.json configuration file with following sections:

### Method Parameters

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

### Dust
- kappa &mdash; The extinction opacity
- albedo &mdash; Single scattering albedo
- hgg &mdash; Henyey-Greenstein phase function anisotropy
- pl &mdash; Peak linear polarization
- pc &mdash; Peak circular polarization
- sc &mdash; Skew Factor

### Geometry

The geometry of the scattering matter is described as a set of geometric shapes. 
On the top level, there are should be flared disk, sphere envelope or max/sum list.
All angles of this section are measured in degrees, all distances are supposed to be in astronomy units.

#### Flared Disk

The density of the flared disk is described by the following formulae:

![flared disk](./docs_src/images/flared_disk.svg)

where the scale height

![flared disk height](./docs_src/images/flared_disk_h.svg)

and the radial coordinate in the disk plane

![flared disk radius](./docs_src/images/flared_disk_r.svg).

The model parameters are:
- rInner &mdash; the inner disk radius
- rOuter &mdash; the outer disk radius
- rho0 &mdash; the density at the disk midplane at a radius r0
- h0 &mdash; the disk scale height at a radius r0
- r0 &mdash; the radius, where h0 and rho0 are defined
- alpha &mdash; the radial density exponent
- beta &mdash; the flaring power

##### Safier Wind

One can add the disk wind to the flared disk.
The density of the disk with the wind is the maximum of the wind and disk densities.
The density of the Safier wind is
 
![safier wind](./docs_src/images/safier_wind.svg)

where _&chi; = z / r_ is the dimensionless height above the disk plane and _&rho;<sub>0</sub>_ is the wind density on
the disk surface at distance 1 AU from the star.

![safier wind](./docs_src/images/safier_wind_rho0.svg)

The function _&eta; (&chi;)_ can be obtained by solving the gas-dynamic equations.
_&xi;'<sub>0</sub>_ and _&psi;<sub>0</sub>_ are wind model parameters, which are defined in Safier papers for a list of models.

The parameters of the wind are:
- model &mdash; the model of the wind from Safier papers. Should be B, C, D, I, E, F or G
- mOut &mdash; the mass outflow rate in solar masses per year
- mStar &mdash; the stellar mass in solar masses
- h0 &mdash; the dimensionless height from which the wind begins
- rMin &mdash; the inner radius of the wind formation region
- rMax &mdash; the outer radius of the wind formation region

#### Sphere Envelope

The density of the sphere envelope is described by the equation

![sphere envelope](./docs_src/images/sphere_envelope.svg)

where the radius is
    
![sphere envelope radius](./docs_src/images/sphere_envelope_r.svg).

The model parameters are:
- rInner &mdash; the inner sphere radius
- rOuter &mdash; the outer sphere radius
- rho0 &mdash; the density at a radius r0
- r0 &mdash; the radius, where rho0 is defined
- alpha &mdash; the radial density exponent
      
#### Sum / Max

The list of other geometry shapes (including other lists).
The density of the sum list is the sum of densities of all elements.
The density of the max list is the maximum density of all elements.

#### Translation

The Flared Disk and the Sphere Envelope may contain additional Translation section.
It describes a rotation and a translation of the shape in the global coordinate system.
We use the Euler angles for rotation.
So there are six parameters: intrinsicRotation, nutation, precession, x, y, z.
All transformations are applied in the global coordinate system, thus the order of operations is:
1. Intrinsic Rotation (rotation around _z_-axis)
2. Nutation (rotation around _x_-axis)
3. Precession (rotation around _z_-axis)
4. (x, y, z) vector translation

Every value may be omitted (it will be treated as zero).

### Stars
Is a list of sources with 4 parameters:
- x, y, z &mdash; source coordinates
- l &mdash; source luminosity

```yaml
"stars": [
    {
      "x": 0.0,
      "y": 0.0,
      "z": 0.0,
      "l": 1.0
    }
  ]
```

### Observers

We use spherical coordinates [_&theta;_, _&phi;_] to describe the direction towards each observer.
_&theta;_ is a zenith distance (an angle between positive _z_-axis direction and observers direction), measured from 0 to 180&deg;.
_&phi;_ is an azimuth, measured counterclockwise in _xy_ plane from the positive _x_-axis direction.
The azimuth may be in the range from -180&deg; to 180&deg; or from 0 to 360&deg;.

The third parameter of the observer is _rimage_.
It is a radius of the visible area in astronomy units.
Both axes of the image will be from -_rimage_ to _rimage_.
We specify _rimage_ once for all observers.

There are three ways to specify observer directions:

#### Manual

One should specify the list of coordinates _&phi;_ and _&theta;_ for each observer position.

#### Parallel

All observers are evenly distributed on the circle with a constant _&theta;_.
In this configuration, one should specify the number of observers _numberOfObservers_ and _&theta;_.
The first observer has _&phi;_ equal to zero.

#### Meridian

All observers are evenly distributed on the meridian with a constant _&phi;_.
In this configuration, one should specify the number of observers _numberOfObservers_ and _&phi;_.
If _numberOfObservers_ is the odd one of the observers has _&theta;_ equal to 90&deg;.
This configuration does not include poles.

#### Example

You can combine all three ways in one problem.
The JSON configuration may be specified in the following way:

```yaml
"observers": {
    "rimage": 800.0,
    "manual": [
      {
        "phi": 45.0,
        "theta": 45.0
      },
      {
        "phi": 60.0,
        "theta": 120.0
      }
    ],
    "parallel": {
      "numberOfObservers" : 10,
      "theta" : 90.0
    },
    "meridian": {
      "numberOfObservers" : 2,
      "phi" : 0.0
    }
  }
```

## Used third-party libraries

* **[nlohmann/json](https://github.com/nlohmann/json)** to parse a configuration file
* **[catch2](https://github.com/catchorg/Catch2)** for unit testing

