#include "FlaredDisk.hpp"
#include <cmath>

double FlaredDisk::density(double const x, double const y, double const z) const
{
    double const r2 = x*x + y*y;
    double const r  = std::sqrt(r2);

    // Disk Geometry
    if(( r >= rInner_ ) && ( r <= rOuter_ ))
    {
        double const h = h0_ * std::pow(r/r0_, beta_);
        return rho0_ * std::pow(r0_ / r, alpha_) * std::exp(-0.5*z*z / (h*h)); // rho in g/cm^3
    } else {
        return 0.0;
    }
}