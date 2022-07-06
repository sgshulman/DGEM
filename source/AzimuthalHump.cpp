#include "AzimuthalHump.hpp"
#include "DebugUtils.hpp"
#include <cmath>

AzimuthalHump::AzimuthalHump(
    double const h,
    double const r,
    double const sigma2,
    double const sigma2azimuthalForward,
    double const sigma2azimuthalBackward)
    : h_{ h }
    , r_{ r }
    , sigma2_{ sigma2 }
    , sigma2azForward_{ sigma2azimuthalForward }
    , sigma2azBackward_{ sigma2azimuthalBackward > 0.0 ? sigma2azimuthalBackward : sigma2azimuthalForward }
{
    DATA_ASSERT(h_ >= 0., "h (the azimuthal hump height) must be non negative.");
    DATA_ASSERT(r_ > 0., "r (the azimuthal hump distance from the disk axis) must be positive.");
    DATA_ASSERT(sigma2_ > 0., "sigma2 for the azimuthal hump must be positive.");
    DATA_ASSERT(sigma2azForward_ > 0., "sigma2azimuthal for the azimuthal hump must be positive.");
}


double AzimuthalHump::hump(double value, const Vector3d &position) const
{
    double const r = std::sqrt(position.x() * position.x() + position.y() * position.y());
    double const az = std::atan2(position.y(), position.x());
    double const azExpBase = az > 0 ? (az * az / (2*sigma2azForward_)) : (az * az / (2*sigma2azBackward_));

    return value * (1.0 + h_ *
        std::exp( - (r-r_)*(r-r_) / (2*sigma2_) )*
        std::exp( - azExpBase ));
}
