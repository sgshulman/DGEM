#include "FractalCloud.hpp"
#include "DebugUtils.hpp"
#include "LEcuyer.hpp"
#include <vector>

namespace
{
    std::vector<Vector3d> buildSmallCubes(
        double const d,
        std::uint32_t const cubeDotNumber,
        std::uint32_t const oldDotNumber,
        std::vector<Vector3d> const& dots,
        LEcuyer* ran)
    {
        std::vector<Vector3d> newDots(cubeDotNumber * oldDotNumber);

        for (std::uint64_t i = 0; i != oldDotNumber; ++i)
        {
            for (std::uint64_t j = 0; j != cubeDotNumber; ++j)
            {
                double const z = dots.at(i).z() + ran->Get()*2*d - d;
                double const y = dots.at(i).y() + ran->Get()*2*d - d;
                double const x = dots.at(i).x() + ran->Get()*2*d - d;
                newDots.at(i*cubeDotNumber+j) = Vector3d(x, y, z);
            }
        }

        return newDots;
    }
}


FractalCloud::FractalCloud(
    std::uint32_t const n,
    double const max,
    double const dCube,
    double const rho0,
    std::uint32_t const dotsN,
    std::int32_t const seed)
        : n_(n)
        , max_(max)
{
    DATA_ASSERT(rho0 > 0, "rho0 (density of fractal cloud) must be positive.");
    DATA_ASSERT(dCube > 0, "dCube (fractal dimension) must be positive.");

    const double l = 2 * max;
    LEcuyer ran(seed);

    std::vector<Vector3d> dots0;
    dots0.reserve(dotsN);

    for (std::uint64_t i = 0; i != dotsN; ++i)
    {
        double const z = ran.Get()*l - max_;
        double const y = ran.Get()*l - max_;
        double const x = ran.Get()*l - max_;
        dots0.emplace_back(x, y, z);
    }

    const double delta = std::pow(1.0*dotsN, 1.0/dCube);

    const double d1 = l/(2*delta);
    std::vector<Vector3d> dots1 = buildSmallCubes(d1, dotsN, dotsN, dots0, &ran);

    const double d2 = d1/(2*delta);
    const std::uint32_t dotsN2 = dotsN * dotsN;
    std::vector<Vector3d> dots2 = buildSmallCubes(d2, dotsN, dotsN2, dots1, &ran);

    const double d3 = d2/(2*delta);
    const std::uint32_t dotsN3 = dotsN2 * dotsN;
    std::vector<Vector3d> dots3 = buildSmallCubes(d3, dotsN, dotsN3, dots2, &ran);

    density_ = new double[n_*n_*n_]();
    const std::uint32_t dotsN4 = dotsN3 * dotsN;

    for (std::uint64_t i = 0; i != dotsN4; ++i)
    {
        double x = dots3[i].x();
        while (x < -max_) x += l;
        while (x >  max_) x -= l;

        double y = dots3[i].y();
        while (y < -max_) y += l;
        while (y >  max_) y -= l;

        double z = dots3[i].z();
        while (z < -max_) z += l;
        while (z >  max_) z -= l;

        auto const cntx = std::uint32_t((x + max_)/(2*max_)*n_);
        auto const cnty = std::uint32_t((y + max_)/(2*max_)*n_);
        auto const cntz = std::uint32_t((z + max_)/(2*max_)*n_);

        density_[cntx+cnty*n_+cntz*n_*n_] += rho0;
    }
}


FractalCloud::~FractalCloud()
{
    delete[] density_;
}


double FractalCloud::density(Vector3d const& position) const
{
    auto const cntx = std::uint32_t((position.x() + max_)/(2*max_)*n_);
    auto const cnty = std::uint32_t((position.y() + max_)/(2*max_)*n_);
    auto const cntz = std::uint32_t((position.z() + max_)/(2*max_)*n_);
    return density_[cntx+cnty*n_+cntz*n_*n_];
}
