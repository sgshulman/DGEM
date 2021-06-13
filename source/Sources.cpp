#include "Sources.hpp"

#include <iostream>
#include "Observer.hpp"
#include "CartesianGrid.hpp"
#include "IRandomGenerator.hpp"
#include "MathUtils.hpp"

namespace
{
    inline std::pair<Vector3d, std::uint64_t> starSurface(
        SphereSource const& source,
        Direction3d const& direction,
        IGridCRef grid)
    {
        Photon innerPh(source.pos(), source.cellId(), direction,1.0,0);
        grid->movePhotonAtDistance(innerPh, source.radius());
        return { innerPh.pos(), innerPh.cellId() };
    }

    inline Direction3d randomDirection(IRandomGenerator* ran)
    {
        double const v = ran->Get();
        double const u = ran->Get();
        return { 2.0 * PI * u, std::acos(2 * v - 1.0) };
    }
}


Photon Sources::emitRandomPhoton(IGridCRef grid, IRandomGenerator* ran)
{
    if (currentSource_ == pointSources_.size() + sphereSources_.size())
    {
        return {Vector3d{}, 0, Direction3d{}, 0.0, std::numeric_limits<std::uint32_t>::max()};
    }

    std::uint32_t const sourceId{ currentSource_ };
    std::uint64_t const photonId{ photonId_ };
    ++photonId_;

    if(photonId % 10000 == 0)
    {
        std::cout << "Sources: " << sourceId + 1 << "/" << pointSources_.size() + sphereSources_.size() <<
                  ". Photons: " << photonId << "/" << photonsNumber_ << std::endl;
    }

    if (photonId_ == photonsNumber_)
    {
        ++currentSource_;
        photonId_ = 0;

        double const nextLuminosity = currentSource_ < pointSources_.size()
            ? pointSources_[currentSource_].luminosity()
            : currentSource_ - pointSources_.size() < sphereSources_.size()
            ? sphereSources_[currentSource_ - pointSources_.size()].luminosity() : 0;

        photonsNumber_ = (std::uint64_t) (parameters_.num_photons_ * nextLuminosity / totlum_);
    }

    if (sourceId < pointSources_.size())
    {
        return {
            pointSources_[sourceId].pos(),
            pointSources_[sourceId].cellId(),
            randomDirection(ran),
            1.0,
            1};
    }

    std::uint32_t const sphereSourceId = sourceId - static_cast<std::uint32_t>(pointSources_.size());

    // compute point location
    Direction3d const positionDirection{ randomDirection(ran) };
    auto const position = starSurface(sphereSources_[sphereSourceId], positionDirection, grid);

    double const v = ran->Get();
    double const u = ran->Get();
    Direction3d localDirection{ 2.0 * PI * u, std::sqrt(v) };

    return { position.first, position.second, localDirection.rotate(positionDirection), 1.0, 1};
}


Photon Sources::emitDgemPhoton(IGridCRef grid)
{
    if (currentSource_ == pointSources_.size() + sphereSources_.size())
    {
        return {Vector3d{}, 0, Direction3d{}, 0.0, std::numeric_limits<std::uint32_t>::max()};
    }

    std::uint32_t const sourceId{ currentSource_ };
    std::uint64_t const photonId{ photonId_ };

    ++photonId_;

    if(photonId % 1000 == 0)
    {
        std::cout << "Sources: " << sourceId + 1 << "/" << pointSources_.size() + sphereSources_.size() <<
                  ". Photons: " << photonId << "/" << photonsNumber_ << std::endl;
    }

    if (photonId_ == photonsNumber_)
    {
        ++currentSource_;
        photonsNumber_ = primaryDir_.number();
        photonId_ = 0;
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
        auto const position = starSurface(sphereSources_[sphereSourceId], sphereDir_.direction(pointId_), grid);
        pointPosition_ = position.first;
        pointCellId_ = position.second;
        ++pointId_;
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
            if ((*observers)[io].inFov(pointSources_[is].pos()))
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
    }

    // Direct photon loop. Loop over sphere sources and over points on them,
    // weight photons by W=ph*exp(-tau1)/2pi/number_of_directions
    for (std::uint64_t is=0; is!=sphereSources_.size(); ++is)
    {
        auto const nph = std::uint64_t(parameters_.num_photons_ * sphereSources_[is].luminosity() / totlum_);

        for (std::uint64_t ip=0; ip!=sphereDir_.number(); ++ip)
        {
            auto const position = starSurface(sphereSources_[is], sphereDir_.direction(ip), grid);

            for (std::uint64_t io = 0; io != observers->size(); ++io)
            {
                double cosTheta = sphereDir_.direction(ip) * (*observers)[io].direction().vector();

                if (cosTheta >= 0 && (*observers)[io].inFov(position.first))
                {
                    // Set photon location, grid cell, and direction of observation
                    Photon ph(position.first, position.second, (*observers)[io].direction(), 1.0, 0);

                    // Find optical depth, tau1, to edge of grid along viewing direction
                    double tau1 = grid->findRealOpticalDepth(sphereSources_[is].pos(),
                                                             (*observers)[io].direction().vector());

                    // direct photon weight is exp(-tau1)/2pi/number_of_directions
                    ph.weight() = nph * std::exp(-tau1) / 2.0 / PI / sphereDir_.number() * cosTheta;

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

    for (std::uint64_t is=0; is!=sphereSources_.size(); ++is)
    {
        std::cout << "Optical depths from the sphereSource center #" << is << "." << std::endl;

        for (std::uint64_t io=0; io!=observers->size(); ++io)
        {
            Observer const& obs = (*observers)[io];

            // Set photon location, grid cell, and direction of observation
            Photon const ph(sphereSources_[is].pos(), sphereSources_[is].cellId(), obs.direction(), 1.0, 0);

            // Find optical depths to edge of the grid along the line of sight
            double const tau = grid->findOpticalDepth(ph);
            double const realTau = grid->findRealOpticalDepth(sphereSources_[is].pos(), obs.direction().vector());

            std::cout << "Observer #" << io << " (" << degrees(obs.phi()) << ", "
                      << degrees(obs.theta()) << ")\t" << tau << "\t" << realTau << std::endl;
        }
    }
}
