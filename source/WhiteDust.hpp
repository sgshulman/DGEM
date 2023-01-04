#ifndef WHITE_DUST_HPP_
#define WHITE_DUST_HPP_

#include "IDust.hpp"

// The dust with Henyey and Greenstein (1941) phase function and
// White (1979) approximation for polarization functions.
class WhiteDust final: public IDust
{
public:
    WhiteDust(
        double const albedo,
        double const hgg,
        double const pl,
        double const pc,
        double const sc)
    : albedo_{ albedo }
    , hgg_{ hgg }
    , hgg2_{ hgg_ * hgg_ }
    , pl_{ pl }
    , pc_{ pc }
    , sc_{ sc }
    {}

    void scatteringMatrixElements(
        double &p1,
        double &p2,
        double &p3,
        double &p4,
        double cosTheta) const override;

    double fraction(double cosTheta) const override;
    double cosRandomTheta(double v) const override;

    double albedo() const override
    {	return albedo_;	}

private:
    double const albedo_;
    double const hgg_;
    double const hgg2_;
    double const pl_;
    double const pc_;
    double const sc_;
};

#endif