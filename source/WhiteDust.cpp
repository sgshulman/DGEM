#include <cmath>
#include "WhiteDust.hpp"
#include "MathUtils.hpp"

void WhiteDust::scatteringMatrixElements(
    double &p1, double &p2, double &p3, double &p4, double cosTheta) const
{
    // Calculate the elements of the phase matrix for a
    // simple representation of the mrn dust mixture using the algorithms
    // for the ultraviolet region due to Richard L. White Ap.J. 229, 954, 1979.
    if(cosTheta >= 1.0) cosTheta = 1.0;
    if(cosTheta <=-1.0) cosTheta =-1.0;

    double const cos2Theta = cosTheta * cosTheta;

    p1 = (1. - hgg2_) / std::pow(1. + hgg2_ - 2. * hgg_ * cosTheta, 1.5);
    p2 = -pl_ * p1 * (1. - cos2Theta)/(1. + cos2Theta);
    p3 = 2. * p1 * cosTheta / (1. + cos2Theta);

    double const theta = std::acos(cosTheta);
    double const f = std::exp(-7.0 * theta / PI);
    double const thetaF = theta * (1. + 3.13 * sc_ * f);
    double const c = std::cos(thetaF);
    double const c2 = c * c;
    p4 = -pc_ * p1 * (1 - c2) / (1 + c2);
}


double WhiteDust::fraction(double const cosTheta) const
{
    return (1. - hgg2_) / std::pow(1. + hgg2_ - 2. * hgg_ * cosTheta, 1.5);
}


double WhiteDust::cosRandomTheta(double const v) const
{
    double const fraction = (1. - hgg2_) / (1. - hgg_ + 2 * hgg_ * v);
    return (1. + hgg2_ - fraction * fraction) / (2. * hgg_);
}
