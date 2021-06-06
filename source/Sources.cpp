#include "Sources.hpp"

#include <iostream>
#include "Observer.hpp"
#include "CartesianGrid.hpp"
#include "MathUtils.hpp"
#include "Random.hpp"

Photon Sources::emitPhoton(IGridCRef grid, Random* ran)
{
    if (currentSource_ == pointSources_.size() + sphereSources_.size())
    {
        return {Vector3d{}, 0, Direction3d{}, 0.0, std::numeric_limits<std::uint32_t>::max()};
    }

    std::uint32_t const sourceId{ currentSource_ };
    std::uint64_t const photonId{ photonId_ };

    ++photonId_;

    std::uint64_t const period = parameters_.useMonteCarlo_ ? 10000 : 1000;

    if(photonId % period == 0)
    {
        std::cout << "Sources: " << sourceId + 1 << "/" << pointSources_.size() + sphereSources_.size() <<
                  ". Photons: " << photonId << "/" << photonsNumber_ << std::endl;
    }

    if (photonId_ == photonsNumber_)
    {
        ++currentSource_;

        double const nextLuminosity = currentSource_ < pointSources_.size()
            ? pointSources_[currentSource_].luminosity()
            : sphereSources_[currentSource_ - pointSources_.size()].luminosity();

        photonsNumber_ = parameters_.useMonteCarlo_
                        ? (std::uint64_t) (parameters_.num_photons_ * nextLuminosity / totlum_)
                        : primaryDir_.number();

        photonId_ = 0;
    }

    if (parameters_.useMonteCarlo_)
    {
        if (sourceId < pointSources_.size())
        {
            double const v = ran->Get();
            double const u = ran->Get();

            return {
                pointSources_[sourceId].pos(),
                pointSources_[sourceId].cellId(),
                Direction3d(2.0 * PI * u, std::acos(2 * v - 1.0)),
                1.0,
                1};
        }

        std::uint32_t const sphereSourceId = sourceId - static_cast<std::uint32_t>(pointSources_.size());

        // compute point location
        double const v = ran->Get();
        double const u = ran->Get();

        Photon innerPh(
            sphereSources_[sphereSourceId].pos(),
            sphereSources_[sphereSourceId].cellId(),
            Direction3d(2.0 * PI * u, std::acos(2 * v - 1.0)),
            1.0,
            0);

        grid->movePhotonAtDistance(innerPh, sphereSources_[sphereSourceId].radius());

        double const v1 = ran->Get();
        double const u1 = ran->Get();

        // TODO: add check for half-sphere
        return {
            innerPh.pos(),
            innerPh.cellId(),
            Direction3d(2.0 * PI * u1, std::acos(2 * v1 - 1.0)),
            1.0,
            1};
    }

    if (sourceId < pointSources_.size())
    {
        return {
            pointSources_[sourceId].pos(),
            pointSources_[sourceId].cellId(),
            primaryDir_.direction(photonId),
            primaryDir_.w(photonId) * pointSources_[sourceId].luminosity() / totlum_,
            1};
    }

    std::uint32_t const sphereSourceId = sourceId - static_cast<std::uint32_t>(pointSources_.size());

    // TODO: add real logic
    if (photonId % sphereSources_.size() == 0)
    {
        Photon innerPh(
            sphereSources_[sphereSourceId].pos(),
            sphereSources_[sphereSourceId].cellId(),
            primaryDir_.direction(photonId),
            1.0,
            0);

        grid->movePhotonAtDistance(innerPh, sphereSources_[sphereSourceId].radius());

        pointPosition_ = innerPh.pos();
        pointCellId_ = innerPh.cellId();
    }

    // TODO: add check for half-sphere
    return {
        pointPosition_,
        pointCellId_,
        primaryDir_.direction(photonId),
        primaryDir_.w(photonId) * sphereSources_[sphereSourceId].luminosity() / totlum_,
        1};
}


void Sources::directPhotons(IGridCRef grid, std::vector<Observer>* observers)
{
    // Direct photon loop.  Loop over sources and weight photons by
    // W=ph*exp(-tau1)/4pi
    for (std::uint64_t is=0; is!=pointSources_.size(); ++is)
    {
        auto const nph = std::uint64_t(parameters_.num_photons_ * pointSources_[is].luminosity() / totlum_);

        for (std::uint64_t io=0; io!=observers->size(); ++io)
        {
            // Set photon location, grid cell, and direction of observation
            Photon ph(pointSources_[is].pos(), pointSources_[is].cellId(), (*observers)[io].direction(), 1.0, 0);

            // Find optical depth, tau1, to edge of grid along viewing direction
            double tau1 = grid->findRealOpticalDepth(pointSources_[is].pos(), (*observers)[io].direction().vector());

            // direct photon weight is exp(-tau1)/4pi
            ph.weight() = nph * std::exp(-tau1) / 4.0 / PI;

            // bin the photon into the image according to its position and
            // direction of travel.
            (*observers)[io].bin(ph);
        }
    }

    // Direct photon loop. Loop over sphere sources and over points on them,
    // weight photons by W=ph*exp(-tau1)/2pi/number_of_directions
    for (std::uint64_t is=0; is!=sphereSources_.size(); ++is)
    {
        auto const nph = std::uint64_t(parameters_.num_photons_ * pointSources_[is].luminosity() / totlum_);

        for (std::uint64_t ip=0; ip!=sphereDir_.number(); ++ip)
        {
            // compute point location
            Photon innerPh(sphereSources_[is].pos(), sphereSources_[is].cellId(), sphereDir_.direction(ip), 1.0, 0);
            grid->movePhotonAtDistance(innerPh, sphereSources_[is].radius());

            for (std::uint64_t io = 0; io != observers->size(); ++io)
            {
                if (sphereDir_.direction(ip) * (*observers)[io].direction().vector() >= 0)
                {
                    // Set photon location, grid cell, and direction of observation
                    Photon ph(innerPh.pos(), innerPh.cellId(), (*observers)[io].direction(), 1.0, 0);

                    // Find optical depth, tau1, to edge of grid along viewing direction
                    double tau1 = grid->findRealOpticalDepth(pointSources_[is].pos(),
                                                             (*observers)[io].direction().vector());

                    // direct photon weight is exp(-tau1)/2pi/number_of_directions
                    ph.weight() = nph * std::exp(-tau1) / 2.0 / PI / sphereDir_.number();

                    // bin the photon into the image according to its position and
                    // direction of travel.
                    (*observers)[io].bin(ph);
                }
            }
        }
    }
}


void Sources::writeObserversOpticalDepths(IGridCRef grid, std::vector<Observer>* observers)
{
    for (std::uint64_t is=0; is!=pointSources_.size(); ++is)
    {
        std::cout << "Optical depths from the pointSource #" << is << "." << std::endl;

        for (std::uint64_t io=0; io!=observers->size(); ++io)
        {
            Observer const& obs = (*observers)[io];

            // Set photon location, grid cell, and direction of observation
            Photon const ph(pointSources_[is].pos(), pointSources_[is].cellId(), obs.direction(), 1.0, 0);

            // Find optical depths to edge of the grid along the line of sight
            double const tau = grid->findOpticalDepth(ph);
            double const realTau = grid->findRealOpticalDepth(pointSources_[is].pos(), obs.direction().vector());

            std::cout << "Observer #" << io << " (" << degrees(obs.phi()) << ", "
                << degrees(obs.theta()) << ")\t" << tau << "\t" << realTau << std::endl;
        }
    }
}
