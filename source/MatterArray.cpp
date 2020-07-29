#include "MatterArray.hpp"

double MatterArray::density(Vector3d const& position) const
{
    double rho = 0.;

    for (const auto& matter : matterArray_)
    {
        rho = (*accumulator_)(rho, matter->density(position));
    }

    return rho;
}
