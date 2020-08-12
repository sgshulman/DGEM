#ifndef CARTESIAN_GRID_HPP_
#define CARTESIAN_GRID_HPP_

#include "IGrid.hpp"
#include "Predefines.hpp"
#include <cstdint>
#include <memory>

// Cartesian grid
class CartesianGrid : public IGrid
{
    public:
        CartesianGrid(double xmax, double ymax, double zmax, double kappa,
                      std::uint32_t nx, std::uint32_t ny, std::uint32_t nz, IMatterCPtr matter);

        ~CartesianGrid() override
        {
            delete[] rhokappa_;
        }

        double calculateRealTau(const Vector3d& v) const override;
        double findOpticalDepth(Photon ph) const override;
        int movePhotonAtDepth(Photon& ph, double tau, double tauold) const override;
        int movePhotonAtRandomDepth(Photon& ph, Random *ran) const override;
        void peeloff(Photon ph, Observer& observer, DustCRef dust) const override;
        double computeMatterMass() const override;
        std::uint64_t cellId(const Vector3d& position) const override;

    private:
        double maxDistance(Photon const& ph) const;
        std::pair<double, std::uint64_t> cellDistance(Photon& ph, double delta, Vector3d const& phDirInv) const;

        std::uint32_t nx_, ny_, nz_;
        double *rhokappa_{ nullptr };
        double xmax_, ymax_, zmax_;
        double const xCellSize_, yCellSize_, zCellSize_;
        double const xCellSizeInv_, yCellSizeInv_, zCellSizeInv_;
        double minrho_;
        double const kappa_;
        IMatterCPtr matter_;
};

#endif
