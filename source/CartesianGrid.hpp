#ifndef CARTESIAN_GRID_HPP_
#define CARTESIAN_GRID_HPP_

#include "IGrid.hpp"
#include <cstdint>
#include <memory>

// Cartesian grid
class CartesianGrid : public IGrid
{
    public:
        CartesianGrid(double xmax, double ymax, double zmax, double kappa,
                      uint32_t nx, uint32_t ny, uint32_t nz, IMatterCPtr matter);

        ~CartesianGrid() override
        {
            delete[] rhokappa_;
        }

        double findOpticalDepth(Photon ph) const override;
        int movePhotonAtDepth(Photon& ph, double tau, double tauold) const override;
        int movePhotonAtRandomDepth(Photon& ph, Random *ran) const override;
        void peeloff(Photon ph, Observer& obs, DustCRef dust) const override;
        double computeMatterMass() const override;

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
