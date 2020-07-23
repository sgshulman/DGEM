#ifndef DIRECTION_3D_HPP_
#define DIRECTION_3D_HPP_

#include "Vector3d.hpp"
#include "MathUtils.hpp"

class Direction3d
{
public:
    Direction3d() = default;

    Direction3d(Vector3d const vector)
        : vector_{ vector.normalized() }
        , phi_{ std::atan2(vector_.y(), vector_.x()) }
        , sinTheta_{ std::sqrt(vector_.x()*vector_.x() + vector_.y()*vector_.y()) }
    {}

    Direction3d(double const phi, double const theta)
        : phi_{ normAngle(phi) }
        , sinTheta_{ std::sin(theta) }
    {
        vector_ = {sinTheta_ * std::cos(phi_), sinTheta_ * std::sin(phi_), std::cos(theta) };
    }

    inline Vector3d vector() const
    {   return vector_; }
    inline double x() const
    { return vector_.x(); }
    inline double y() const
    { return vector_.y(); }
    inline double z() const
    { return vector_.z(); }
    inline double phi() const
    {	return phi_;	}
    inline double cosTheta() const
    {	return vector_.z();	}
    inline double sinTheta() const
    {	return sinTheta_;	}

private:
    Vector3d vector_;
    double phi_{ 0. };
    double sinTheta_{ 0. };
};

#endif //DIRECTION_3D_HPP_
