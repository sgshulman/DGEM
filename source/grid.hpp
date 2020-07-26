#ifndef GRID_HPP_
#define GRID_HPP_

#include <cstdint>
#include <memory>

#include "Predefines.hpp"

// Cartesian grid
class Grid
{
    public:
        Grid(double xmax, double ymax, double zmax, double kappa,
            uint32_t nx, uint32_t ny, uint32_t nz, IMatterCPtr matter);

        ~Grid()
        {
            delete[] rhokappa_;
        }

        Grid(Grid const&) = delete;
        Grid& operator=(Grid const&) = delete;

        double findOpticalDepth(Photon ph, double delta=-0.001) const;
        int movePhotonAtDepth(Photon& ph, double tau, double tauold=0.0, double delta=-0.001) const;
        int movePhotonAtRandomDepth(Photon& ph, Random *ran, double delta=-0.001) const;
        void peeloff(Photon ph, Observer& obs, DustCRef dust) const;

    private:
        double maxDistance(Photon const& ph) const;
        double cellDistance(Photon& ph, double delta) const;

        uint32_t nx_, ny_, nz_;
        double *rhokappa_{ nullptr };
        double xmax_, ymax_, zmax_;
        double minrho_;
        IMatterCPtr matter_;
};

#endif
