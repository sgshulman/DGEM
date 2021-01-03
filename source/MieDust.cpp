#include <cmath>
#include "MieDust.hpp"

void MieDust::scatteringMatrixElements(
    double &p1, double &p2, double &p3, double &p4, double cosTheta) const
{
    (void)p1;
    (void)p2;
    (void)p3;
    (void)p4;
    (void)cosTheta;
}


double MieDust::fraction(double const cosTheta) const
{
    (void)cosTheta;
    return 0.0;
}


double MieDust::cosRandomTheta(double const v) const
{
    (void)v;
    return 0.0;
}
