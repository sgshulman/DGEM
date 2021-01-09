#include "KurosawaWind.hpp"
#include "DebugUtils.hpp"
#include <cmath>

KurosawaWind::KurosawaWind(
        double const d,
        double const p,
        double const mOut,
        double const mStar,
        double const rInner,
        double const rOuter,
        double const rScale,
        double const windAccelerationRate,
        double const terminalV,
        double const soundSpeed)
    : d_{ d }
    , mStar_{ mStar }
    , rScale_{ rScale }
    , windAccelerationRate_{ windAccelerationRate }
    , terminalV_{ terminalV }
    , soundSpeed_{ soundSpeed }
{
    DATA_ASSERT(rInner > 0., "rInner (inner radius of the wind formation region) must be positive.");
    DATA_ASSERT(rOuter > 0., "rOuter (outer radius of the wind formation region) must be positive.");
    DATA_ASSERT(rOuter > rInner, "rOuter (outer radius) must be greater than rInner (inner radius of the wind formation region).");

    cosInner_ = d / std::sqrt(d * d + rInner * rInner);
    cosOuter_ = d / std::sqrt(d * d + rOuter * rOuter);
    mOutNormalization_ = (p + 1.) * mOut / (std::pow(rOuter, p + 1.) - std::pow(rInner, p + 1.));
}


double KurosawaWind::density(Vector3d const& position) const
{
    (void)position;
    return 0.0;
}
