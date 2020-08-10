#ifndef SOURCES_HPP_
#define SOURCES_HPP_

#include <cstdint>
#include "Vector3d.hpp"
#include "Photon.hpp"
#include "Predefines.hpp"


// one point source of photons
class PointSource
{
public:
    PointSource()
        : cellId_{}
        , lum_{}
    {}

    PointSource(Vector3d const& pos, std::uint32_t cellId, double lum)
        : pos_(pos)
        , cellId_{ cellId }
        , lum_(lum)
    {};

    PointSource(PointSource const& other) = delete;
    PointSource& operator=(PointSource const& other) = delete;

    PointSource(PointSource&& other) = default;
    PointSource& operator=(PointSource&& other) = default;

    Vector3d const& pos() const
    {	return pos_; }

    double luminosity() const
    {	return lum_; }

    uint32_t cellId() const
    {   return cellId_; }

private:
    Vector3d	pos_;
    std::uint32_t    cellId_;
    double		lum_;
};

struct SourceParameters
{
    bool useMonteCarlo_;
    uint64_t num_photons_;
    uint32_t PrimaryDirectionsLevel_;
};

// all sources of photons
class Sources
{
public:
    Sources(SourceParameters parameters, std::vector<PointSource> pointSources)
        : parameters_{ parameters}
        , pointSources_{ std::move(pointSources) }
        , totlum_{ 0. }
        , primaryDir_{ parameters_.useMonteCarlo_ ? 1 : parameters_.PrimaryDirectionsLevel_ }
        , currentSource_{ 0 }
        , photonId_{ 0 }
    {
        if (!parameters.useMonteCarlo_)
        {
            parameters_.num_photons_ = primaryDir_.number();
        }

        for (const auto& pointSource : pointSources_)
        {
            totlum_ += pointSource.luminosity();
        }

        photonsNumber_ = parameters_.useMonteCarlo_
                         ? (uint64_t) (parameters_.num_photons_ * pointSources_.at(0).luminosity() / totlum_)
                         : primaryDir_.number();
    }

    ~Sources() = default;

    Sources(Sources const &) = delete;
    Sources& operator=(Sources const&) = delete;

    Photon emitPhoton(Random* ran);
    void directPhotons(IGridCRef grid, std::vector<Observer>* observers);

    size_t num_photons() const
    {   return parameters_.num_photons_; }

private:
    SourceParameters parameters_;
    std::vector<PointSource> const pointSources_;
    double	 totlum_;

    Directions primaryDir_;
    uint32_t currentSource_;
    uint64_t photonId_;
    uint64_t photonsNumber_;
};

#endif
