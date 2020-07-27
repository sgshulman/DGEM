#include "MatterArray.hpp"

double MatterArray::density(double x, double y, double z) const
{
    double rho = 0.;

    for (const auto& matter : matterArray_)
    {
        rho = (*accumulator_)(rho, matter->density(x, y, z));
    }

    return rho;
}
