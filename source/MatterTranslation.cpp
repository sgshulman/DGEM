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

    Vector3d const ir{
        c.x() * cosIntrinsicRotation_ + c.y() * sinIntrinsicRotation_,
        -c.x() * sinIntrinsicRotation_ + c.y() * cosIntrinsicRotation_,
        c.z() };

    Vector3d const n{
        ir.x(),
        ir.y() * cosNutation_ + ir.z() * sinNutation_,
        -ir.y() * sinNutation_ + ir.z() * cosNutation_};

    return {
        n.x() * cosPrecession_ + n.y() * sinPrecession_,
        -n.x() * sinPrecession_ + n.y() * cosPrecession_,
        n.z() };
}
