#include "FlaredDisk.hpp"
#include "DebugUtils.hpp"
#include <algorithm>
#include <cmath>

FlaredDisk::FlaredDisk(
    double const rInner,
    double const rOuter,
    double const rho0,
    double const h0,
    double const r0,
    double const alpha,
    double const beta,
    IMatterCPtr wind)
    : rInner_{ rInner }
    , rOuter_{ rOuter }
    , rho0_{ rho0 }
    , h0_{ h0 }
    , r0_{ r0 }
    , alpha_{ alpha }
    , beta_{ beta }
    , wind_{ std::move(wind) }
{
    DATA_ASSERT(rInner > 0., "rInner (inner radius of the flared disk) must be positive.");
    DATA_ASSERT(rOuter > 0., "rOuter (outer radius of the flared disk) must be positive.");
    DATA_ASSERT(rOuter > rInner, "rOuter (outer radius of the flared disk) must be greater than rInner (inner radius of the flared disk).");
    DATA_ASSERT(rho0 > 0., "rho0 (density of the flared disk) must be positive.");
    DATA_ASSERT(r0 > 0., "r0 (reference radius of the flared disk) must be positive.");
    DATA_ASSERT(h0 > 0., "h0 (the flared disk scale height) must be positive.");
}


double FlaredDisk::density(double const x, double const y, double const z) const
{
    double const r2 = x*x + y*y;
    double const r  = std::sqrt(r2);

    // Disk Geometry
    if(( r >= rInner_ ) && ( r <= rOuter_ ))
    {
        double const h = h0_ * std::pow(r/r0_, beta_);
        double const diskRho = rho0_ * std::pow(r0_ / r, alpha_) * std::exp(-0.5*z*z / (h*h)); // rho in g/cm^3
        double const windRho = wind_ ? wind_->density(x, y, z) : 0.0;
        return std::max(diskRho, windRho);
    }

    return 0.0;
}
