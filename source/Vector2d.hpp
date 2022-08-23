#ifndef VECTOR_2D_HPP_
#define VECTOR_2D_HPP_

#include <cmath>
#include <limits>

class Vector2d
{
    public:
        Vector2d() = default;

        constexpr Vector2d(double const x, double const y)
            : x_(x)
            , y_(y)
        {}

        inline Vector2d normalized() const
        {
            double const l = std::sqrt(x_*x_ + y_*y_);
            return {x_/l, y_/l};
        }

        inline double norm() const
        {
            return std::sqrt(x_*x_ + y_*y_);
        }

        inline double norm2() const
        {
            return x_*x_ + y_*y_;
        }

        // data
        inline double x() const { return x_; }
        inline double y() const { return y_; }
        inline double& x() { return x_; }
        inline double& y() { return y_; }

    private:
        double x_{ 0. };
        double y_{ 0. };
};


inline Vector2d operator+(Vector2d const& left, Vector2d const& right)
{
    return { left.x()+right.x(), left.y()+right.y() };
}

inline Vector2d operator-(Vector2d const& left, Vector2d const& right)
{
    return { left.x()-right.x(), left.y()-right.y() };
}

inline Vector2d operator*(double const left, Vector2d const& right)
{
    return { left*right.x(), left*right.y() };
}

inline Vector2d operator*(Vector2d const& left, double const right)
{
    return right * left;
}

inline Vector2d operator/(Vector2d const& left, double const right)
{
    return (1. / right) * left;
}

inline double operator*(Vector2d const& left, Vector2d const& right)
{
    return left.x()*right.x() + left.y()*right.y();
}

#endif //VECTOR_2D_HPP_
