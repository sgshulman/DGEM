#ifndef I_GRID_HPP_
#define I_GRID_HPP_

#include "Predefines.hpp"

class IGrid
{
    public:
        IGrid() = default;
        virtual ~IGrid() = default;

        IGrid(IGrid const&) = delete;
        IGrid& operator=(IGrid const&) = delete;

        virtual double findOpticalDepth(Photon ph) const = 0;
        virtual int movePhotonAtDepth(Photon& ph, double tau, double tauold) const = 0;
        virtual int movePhotonAtRandomDepth(Photon& ph, Random *ran) const = 0;
        virtual void peeloff(Photon ph, Observer& observer, DustCRef dust) const = 0;
        virtual double computeMatterMass() const = 0;
        virtual std::uint32_t cellId(const Vector3d& position) const = 0;
};

#endif
