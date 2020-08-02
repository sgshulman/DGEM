#include "RoundHump.hpp"
#include "DebugUtils.hpp"
#include <cmath>

RoundHump::RoundHump(double const h, double const r, double const sigma2)
    : h_{ h }
    , r_{ r }
    , sigma2_{ sigma2 }
{
    DATA_ASSERT(h_ >= 0., "h (the round hump height) must be non negative.");
    DATA_ASSERT(r_ > 0., "r (the round hump distance from the disk axis) must be positive.");
    DATA_ASSERT(sigma2_ > 0., "sigma2 for the round hump must be positive.");
}


double RoundHump::hump(double value, const Vector3d &position) const
{
    return value * (1.0 + h_ * std::exp(
        -((position.x()-r_)*(position.x()-r_) + position.y()*position.y())/(2. * sigma2_)));
}
