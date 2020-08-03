#ifndef FRACTAL_CLOUD_HPP_
#define FRACTAL_CLOUD_HPP_

#include "IMatter.hpp"
#include <cstdint>

class FractalCloud : public IMatter
{
    public:
        FractalCloud(uint32_t N, double max, double dCube, double rho0, uint32_t dotsN, int32_t seed);

        ~FractalCloud() override;

        double density(Vector3d const& position) const override;

    private:
        uint32_t const n_;
        double const max_;
        double *density_;
};

#endif
