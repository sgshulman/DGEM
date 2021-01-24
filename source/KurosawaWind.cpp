#include "KurosawaWind.hpp"
#include "DebugUtils.hpp"
#include "MathUtils.hpp"
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
    , p_{ p }
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
    mOutNormalization_ = (p + 2.) * mOut / (std::pow(rOuter, p + 2.) - std::pow(rInner, p + 2.)) / 2. / PI;
}


double KurosawaWind::density(Vector3d const& position) const
{
    Vector3d const rs{0, 0, position.z() > 0. ? -d_ : d_};
    Vector3d const sPos = position - rs;
    double const D = sPos.norm();
    double const cosDelta = -(rs * sPos) / (d_ * D);

    if (cosDelta > cosInner_ || cosDelta < cosOuter_)
    {
        return 0.0;
    }

    double const w = d_ * std::tan(std::acos(cosDelta));
    double const mOut = 0.281565915 * mOutNormalization_ * std::pow(w, p_); // g s-1 cm-2

    double const vCirc = 4.740571715 * std::sqrt(mStar_ / w); // km s-1
    double const vSound = soundSpeed_ * vCirc;
    double const l = D - d_ / cosDelta;
    double const vPoloidal = 1e5 * (vSound + (terminalV_ - vSound) * std::pow(1. - rScale_ / (l + rScale_), windAccelerationRate_)); // sm s-1

    double const right = d_ / (D * cosDelta);
    return mOut / (vPoloidal * std::abs(cosDelta)) * right * right;
}
