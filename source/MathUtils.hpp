#ifndef MATH_UTILS_HPP_
#define MATH_UTILS_HPP_

#include <cmath>

double const PI	= 3.141592653589793238;

inline double degrees(const double radians)
{
    return radians * 180. / PI;
}

inline double radians(const double degrees)
{
    return degrees * PI / 180.;
}

inline double normAngle(const double radians)
{
    double const shift = std::floor(radians / (2. * PI));
    return radians - 2. * PI * shift;
}

inline double clamp(double value, double lower, double high)
{
    if (value < lower)
    {
        return lower;
    }
    if (value > high)
    {
        return high;
    }
    return value;
}

#endif
