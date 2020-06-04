#ifndef VECTOR_3D_HPP_
#define VECTOR_3D_HPP_

#include <cmath>

class Vector3d
{
    public:
        Vector3d()
            : x_(0.0)
            , y_(0.0)
            , z_(0.0)
        {}

        Vector3d(double const x, double const y, double const z)
            : x_(x)
            , y_(y)
            , z_(z)
        {}

        // union vector with spherical coordinates
        Vector3d(double const phi, double const theta)
            : x_(std::cos(phi))
            , y_(std::sin(phi))
            , z_(std::cos(theta))
        {
            const double sinTheta{ std::sin(theta) };
            x_ *= sinTheta;
            y_ *= sinTheta;
        }

        inline Vector3d normalized() const
        {
            double l = std::sqrt(x_*x_ + y_*y_ + z_*z_);
            return {x_/l, y_/l, z_/l};
        }

        inline double norm() const
        {
            return std::sqrt(x_*x_ + y_*y_ + z_*z_);
        }

        // data
        inline double x() const { return x_; }
        inline double y() const { return y_; }
        inline double z() const { return z_; }
        inline double& x() { return x_; }
        inline double& y() { return y_; }
        inline double& z() { return z_; }

    private:
        double x_, y_, z_;
};


inline Vector3d operator+(Vector3d const& left, Vector3d const& right)
{
    return { left.x()+right.x(), left.y()+right.y(), left.z()+right.z() };
}

inline Vector3d operator-(Vector3d const& left, Vector3d const& right)
{
    return { left.x()-right.x(), left.y()-right.y(), left.z()-right.z() };
}

inline Vector3d operator*(double const left, Vector3d const& right)
{
    return { left*right.x(), left*right.y(), left*right.z() };
}

inline Vector3d operator*(Vector3d const& left, double const right)
{
    return right * left;
}

inline Vector3d operator/(Vector3d const& left, double const right)
{
    return (1. / right) * left;
}

inline double operator*(Vector3d const& left, Vector3d const& right)
{
    return left.x()*right.x() + left.y()*right.y() + left.z()*right.z();
}

#endif //VECTOR_3D_HPP_
