#ifndef FRACTAL_CLOUD_HPP_
#define FRACTAL_CLOUD_HPP_

#include "IMatter.hpp"
#include <cstdint>

class FractalCloud final: public IMatter
{
    public:
        FractalCloud(std::uint32_t N, double max, double dCube, double rho0, std::uint32_t dotsN, std::int32_t seed);

        ~FractalCloud() override;

        double density(Vector3d const& position) const override;

    private:
        std::uint32_t const n_;
        double const max_;
        double *density_;
};

#endif
