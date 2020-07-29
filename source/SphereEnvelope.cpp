#include "SphereEnvelope.hpp"
#include "DebugUtils.hpp"
#include "MatterTranslation.hpp"
#include <cmath>

SphereEnvelope::SphereEnvelope(
    double const rInner,
    double const rOuter,
    double const rho0,
    double const r0,
    double const alpha,
    MatterTranslationCPtr translation)
    : rInner_{ rInner }
    , rOuter_{ rOuter }
    , rho0_{ rho0 }
    , r0_{ r0 }
    , alpha_{ alpha }
    , translation_{ std::move(translation) }
{
    DATA_ASSERT(rInner > 0., "rInner (inner radius of the sphere envelope) must be positive.");
    DATA_ASSERT(rOuter > 0., "rOuter (outer radius of the sphere envelope) must be positive.");
    DATA_ASSERT(rOuter > rInner, "rOuter (outer radius) must be greater than rInner (inner radius of the sphere envelope).");
    DATA_ASSERT(rho0 > 0., "rho0 (density of the sphere envelope) must be positive.");
    DATA_ASSERT(r0 > 0., "r0 (reference radius of the sphere envelope) must be positive.");
}


double SphereEnvelope::density(Vector3d const& position) const
{
    Vector3d const pos = translation_ ? (*translation_)(position) : position;
    double const r = pos.norm();

    if (r >= rInner_ && r <= rOuter_)
    {
        return rho0_ * std::pow(r0_ / r, alpha_); // rho in g/cm^3
    }

    return 0.0;
}
