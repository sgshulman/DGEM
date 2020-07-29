#include "MatterTranslation.hpp"
#include <cmath>

MatterTranslation::MatterTranslation(double const precession, double const nutation, double const intrinsicRotation, Vector3d translation)
    : cosPrecession_{ std::cos(precession) }
    , sinPrecession_{ std::sin(precession) }
    , cosNutation_{ std::cos(nutation) }
    , sinNutation_{ std::sin(nutation) }
    , cosIntrinsicRotation_{ std::cos(intrinsicRotation) }
    , sinIntrinsicRotation_{ std::sin(intrinsicRotation) }
    , translation_{ translation }
{}


Vector3d MatterTranslation::operator()(Vector3d const& position) const
{
    Vector3d const c = position - translation_;

    Vector3d const p{
        c.x() * cosPrecession_ + c.y() * sinPrecession_,
        -c.x() * sinPrecession_ + c.y() * cosPrecession_,
        c.z() };

    Vector3d const n{
        p.x(),
        p.y() * cosNutation_ + p.z() * sinNutation_,
        -p.y() * sinNutation_ + p.z() * cosNutation_};

    return {
        n.x() * cosIntrinsicRotation_ + n.y() * sinIntrinsicRotation_,
        -n.x() * sinIntrinsicRotation_ + n.y() * cosIntrinsicRotation_,
        n.z() };
}
