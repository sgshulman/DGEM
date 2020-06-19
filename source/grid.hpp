#ifndef GRID_HPP_
#define GRID_HPP_

#include <cstdint>
#include <memory>

#include "Predefines.hpp"

class FlaredDisk
{
public:
    FlaredDisk(
            double rInner,
            double rOuter,
            double rho0,
            double h0,
            double r0,
            double alpha,
            double beta)
        : rInner_{ rInner }
        , rOuter_{ rOuter }
        , rho0_{ rho0 }
        , h0_{ h0 }
        , r0_{ r0 }
        , alpha_{ alpha }
        , beta_{ beta }
    {}

    double density(double x, double y, double z) const;

private:
    double const rInner_;
    double const rOuter_;
    double const rho0_;
    double const h0_;
    double const r0_;
    double const alpha_;
    double const beta_;
};

// Cartesian grid
class Grid
{
    public:
        Grid(double xmax, double ymax, double zmax, double kappa,
            uint32_t nx, uint32_t ny, uint32_t nz, FlaredDiskCPtr disk);

        ~Grid()
        {
            delete[] rhokappa_;
        }

        Grid(Grid const&) = delete;
        Grid& operator=(Grid const&) = delete;

        double findOpticalDepth(Photon ph, double delta=-0.001) const;
        int movePhotonAtDepth(Photon& ph, double tau, double tauold=0.0, double delta=-0.001) const;
        int movePhotonAtRandomDepth(Photon& ph, double delta=-0.001) const;
        void peeloff(Photon ph, Observer& obs, DustCRef dust) const;

    private:
        double maxDistance(Photon const& ph) const;
        double cellDistance(Photon& ph, double delta) const;

        uint32_t nx_, ny_, nz_;
        double *rhokappa_{ nullptr };
        double xmax_, ymax_, zmax_;
        double minrho_;
        FlaredDiskCPtr disk_;
};

#endif
