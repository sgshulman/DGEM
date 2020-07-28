#include "SphereEnvelope.hpp"
#include "DebugUtils.hpp"
#include <cmath>

SphereEnvelope::SphereEnvelope(
    double const rInner,
    double const rOuter,
    double const rho0,
    double const r0,
    double const alpha)
    : rInner_{ rInner }
    , rOuter_{ rOuter }
    , rho0_{ rho0 }
    , r0_{ r0 }
    , alpha_{ alpha }
{
    DATA_ASSERT(rInner > 0., "rInner (inner radius of the sphere envelope) must be positive.");
    DATA_ASSERT(rOuter > 0., "rOuter (outer radius of the sphere envelope) must be positive.");
    DATA_ASSERT(rOuter > rInner, "rOuter (outer radius) must be greater than rInner (inner radius of the sphere envelope).");
    DATA_ASSERT(rho0 > 0., "rho0 (density of the sphere envelope) must be positive.");
    DATA_ASSERT(r0 > 0., "r0 (reference radius of the sphere envelope) must be positive.");
}


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
