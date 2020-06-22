#include "Sources.hpp"

#include <iostream>
#include "observers.hpp"
#include "grid.hpp"

Photon Sources::emitPhoton()
{
    if (currentSource_ == number_)
    {
        return {Vector3d{}, Direction3d{}, 0.0, std::numeric_limits<uint32_t>::max()};
    }

    uint32_t const sourceId{ currentSource_ };
    uint64_t const photonId{ photonId_ };

    ++photonId_;

    auto const photonsNumber = parameters_.useMonteCarlo_
                               ? (uint64_t) (parameters_.num_photons_ * sources_[sourceId].luminosity() / totlum_)
                               : primaryDir_.number();

    if (photonId_ == photonsNumber)
    {
        ++currentSource_;
        photonId_ = 0;
    }

    if(photonId % 1000 == 0)
    {
        std::cout << "Sources: " << sourceId + 1 << "/" << number_ <<
                  ". Photons: " << photonId << "/" << photonsNumber << std::endl;
    }

    if (parameters_.useMonteCarlo_)
    {
        double const v = ran.Get();
        double const u = ran.Get();

        return {
            sources_[sourceId].pos(),
            Direction3d( 2.0*PI*u, acos( 2*v - 1.0 ) ),
            1.0,
            1};
    }

    return {
        sources_[sourceId].pos(),
        primaryDir_.direction(photonId),
        primaryDir_.w(photonId) * sources_[currentSource_].luminosity() / totlum_,
        1 };
}


void Sources::directPhotons(GridCRef grid, std::vector<Observer>* observers)
{
    // Direct photon loop.  Loop over sources and weight photons by
    // W=ph*exp(-tau1)/4pi
    for (uint64_t is=0; is!=number_; ++is)
    {
        auto const nph = uint64_t(parameters_.num_photons_ * sources_[is].luminosity() / totlum_);

        for (size_t io=0; io!=observers->size(); ++io)
        {
            // Set photon location, grid cell, and direction of observation
            Photon ph(sources_[is].pos(), Direction3d((*observers)[io].pos()), 1.0, 0);

            // Find optical depth, tau1, to edge of grid along viewing direction
            double tau1 = grid->findOpticalDepth(ph);

            // direct photon weight is exp(-tau1)/4pi
            ph.weight() = nph * std::exp(-tau1) / 4.0 / PI;

            // bin the photon into the image according to its position and
            // direction of travel.
            (*observers)[io].Bin(ph);
        }
    }
}
