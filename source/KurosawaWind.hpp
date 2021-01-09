#ifndef KUROSAWA_WIND_HPP_
#define KUROSAWA_WIND_HPP_

#include "IMatter.hpp"

class KurosawaWind : public IMatter
{
public:
    KurosawaWind(
        double d,
        double p,
        double mOut,
        double mStar,
        double rInner,
        double rOuter,
        double rScale,
        double windAccelerationRate,
        double terminalV,
        double soundSpeed);

    ~KurosawaWind() override = default;

    double density(Vector3d const& position) const override;

private:
    double d_;
    double mStar_;
    double cosInner_;
    double cosOuter_;
    double rScale_;
    double windAccelerationRate_;
    double terminalV_;
    double soundSpeed_;
    double mOutNormalization_;
};

#endif
