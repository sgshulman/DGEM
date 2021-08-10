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

        double findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const override;
        double findOpticalDepth(Photon ph) const override;
        double movePhotonAtDistance(Photon& ph, double distance) const override;
        bool movePhotonAtDepth(Photon& ph, double tau, double tauold) const override;
        bool movePhotonAtRandomDepth(Photon& ph, IRandomGenerator *ran) const override;
        void peeloff(Photon ph, Observer& observer, IDustCRef dust) const override;
        void peeloff(Photon ph, Observer &observer, IDustCRef dust, Vector3d const& pos1, Vector3d const& pos2) const override;
        void peeloffHex(Photon ph, Observer& observer, IDustCRef dust, Vector3d const& pos1, Vector3d const& pos2) const override;
        double computeMatterMass() const override;
        double max() const override;
        std::uint64_t cellId(const Vector3d& position) const override;
        bool inside(const Photon& ph) const override;
        void registerSources(SourcesCPtr sources) override;

    private:
        inline bool inside_inner(std::uint64_t cellId) const;

        std::uint32_t const nx_, ny_, nz_;
        std::uint64_t const maxCellId_;
        double *rhokappa_{ nullptr };
        double const xmax_, ymax_, zmax_;
        double const xCellSize_, yCellSize_, zCellSize_;
        double const xCellSizeInv_, yCellSizeInv_, zCellSizeInv_;
        double minrho_;
        double const kappa_;
        IMatterCPtr matter_;
        SourcesCPtr sources_{ nullptr };
};

#endif
