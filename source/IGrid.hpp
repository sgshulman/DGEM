#ifndef I_GRID_HPP_
#define I_GRID_HPP_

#include "Predefines.hpp"
#include <cstdint>

class IGrid
{
    public:
        IGrid() = default;
        virtual ~IGrid() = default;

        IGrid(IGrid const&) = delete;
        IGrid& operator=(IGrid const&) = delete;

        virtual double findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const = 0;
        virtual double findOpticalDepth(Photon ph) const = 0;
        virtual double movePhotonAtDistance(Photon& ph, double distance) const = 0;
        // Returns if photon is inside grid
        virtual bool movePhotonAtDepth(Photon& ph, double tau, double tauold) const = 0;
        virtual bool movePhotonAtRandomDepth(Photon& ph, IRandomGenerator *ran) const = 0;

        virtual void peeloff(Photon ph, Observer& observer, IDustCRef dust) const = 0;
        virtual void peeloff(Photon ph, Observer& observer, IDustCRef dust, Vector3d const& pos1, Vector3d const& pos2) const = 0;
        virtual void peeloffHex(Photon ph, Observer& observer, IDustCRef dust, Vector3d const& pos1, Vector3d const& pos2) const = 0;
        virtual double computeMatterMass() const = 0;
        virtual double max() const = 0;
        virtual std::uint64_t cellId(const Vector3d& position) const = 0;
        virtual bool inside(const Photon& ph) const = 0;
        virtual void registerSources(SourcesCPtr /*sources*/) {};
};

#endif
