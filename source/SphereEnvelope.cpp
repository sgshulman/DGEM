#include "SphereEnvelope.hpp"
#include <cmath>

double SphereEnvelope::density(double const x, double const y, double const z) const
{
    double const r2 = x*x + y*y + z*z;
    double const r = std::sqrt(r2);

    if (r >= rInner_ && r <= rOuter_)
    {
        return rho0_ * std::pow(r0_ / r, alpha_); // rho in g/cm^3
    }

    return 0.0;
}
