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
        virtual int movePhotonAtDepth(Photon& ph, double tau, double tauold) const = 0;
        virtual int movePhotonAtRandomDepth(Photon& ph, Random *ran) const = 0;
        virtual void peeloff(Photon ph, Observer& observer, IDustCRef dust) const = 0;
        virtual double computeMatterMass() const = 0;
        virtual std::uint64_t cellId(const Vector3d& position) const = 0;
};

#endif
