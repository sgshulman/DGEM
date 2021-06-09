#include "Direction3d.hpp"

namespace
{
    inline double sin2cos(double const x)
    {
        if (x * x < 1.)
        {
            return std::sqrt(1. - x * x);
        }

        return 0.;
    }
}


Direction3d Direction3d::rotate(Direction3d const& other) const
{
    // Other.theta is 0 or 180, we rotate around z-axis then result phi = this.phi + other.phi
    if(std::abs(other.sinTheta()) < std::numeric_limits<double>::epsilon())
    {
        if (other.cosTheta() > 0)
        {
            return {phi() + other.phi(), sinTheta(), cosTheta()};
        }
        return {phi() + other.phi(), sinTheta(), -cosTheta()};
    }

    // We consider a spherical triangle with corner angles (A,B,C) and side angles (a,b,c).
    // a is other theta angle (lies on the great circles passing through the z-axis)
    // b is this theta angle
    // c is a result theta angle (lies on the great circles passing through the z-axis)
    // B is (result phi - other phi) (attached to the z-axis)
    // C is this phi angle

    double const cosa = other.cosTheta();
    double const sina = other.sinTheta();

    double const cosb = cosTheta();
    double const sinb = sinTheta();

    double const cosC = cos(phi());
    double sinC;

    if (std::sin(phi()) < 0) // the angle is 2*pi - this phi
    {
        sinC = -sin(phi());
    } else { // the angle is this phi
        sinC = sin(phi());
    }

    bool same_sign;
    double delta;
    if (std::abs(sina) > std::abs(cosa))
    {
        same_sign = sina > 0 && sinb > 0;
        delta = cosb - cosa;
    } else {
        same_sign = cosa > 0 && cosb > 0;
        delta = sinb - sina;
    }

    double cosc;
    double sinc;

    if(same_sign && std::abs(delta) < std::numeric_limits<float>::epsilon() && sinC < std::numeric_limits<float>::epsilon() && cosC > 0)
    {
        if (std::abs(sina) > std::abs(cosa))
        {
            sinc = std::sqrt(delta * delta * (1. + (cosa / sina) * (cosa / sina)) + sina * sinb * sinC * sinC);
        } else {
            sinc = std::sqrt(delta * delta * (1. + (sina / cosa) * (cosa / sina)) + sina * sinb * sinC * sinC);
        }

        cosc = sin2cos(sinc);
    } else {
        cosc = cosa * cosb + sina * sinb * cosC;
        sinc = sin2cos(cosc);
    }

    // result theta is 0 or 180, the result is along the z-axis. phi is undefined
    if(std::abs(sinc) < std::numeric_limits<double>::epsilon())
    {
        if (cosc > 0.)
        {
            return {0., 0.};
        }

        return {0., PI};
    }

    double const cosB = (cosb - cosa * cosc) / (sina * sinc);
    double const sinB = sinC * sinb / sinc;

    // Find final phi values
    if(std::sin(phi()) < 0.)// the top angle is this phi - result phi
    {
        double const finalCosp = cosB * std::cos(other.phi()) + sinB * std::sin(other.phi());
        double const finalSinp = cosB * std::sin(other.phi()) - sinB * std::cos(other.phi());
        return { Vector3d{sinc * finalCosp, sinc * finalSinp, cosc } };
    }
    // the top angle is result phi - this phi
    double const finalCosp = cosB * std::cos(other.phi()) - sinB * std::sin(other.phi());
    double const finalSinp = cosB * std::sin(other.phi()) + sinB * std::cos(other.phi());
    return { Vector3d{sinc * finalCosp, sinc * finalSinp, cosc } };
}
