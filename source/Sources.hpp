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

    PointSource(Vector3d const& pos, std::uint64_t cellId, double lum)
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

    std::uint64_t cellId() const
    {   return cellId_; }

private:
    Vector3d	   pos_;
    std::uint64_t  cellId_;
    double		   lum_;
};

// sphere source of photons
// Limb darkening is ignored
class SphereSource
{
public:
    SphereSource()
        : cellId_{}
        , lum_{}
        , radius_{}
    {}

    SphereSource(Vector3d const& pos, std::uint64_t cellId, double lum, double radius)
        : pos_(pos)
        , cellId_{ cellId }
        , lum_(lum)
        , radius_(radius)
    {};

    SphereSource(SphereSource const& other) = delete;
    SphereSource& operator=(SphereSource const& other) = delete;

    SphereSource(SphereSource&& other) = default;
    SphereSource& operator=(SphereSource&& other) = default;

    Vector3d const& pos() const
    {	return pos_; }

    double luminosity() const
    {	return lum_; }

    double radius() const
    {   return radius_; }

    std::uint64_t cellId() const
    {   return cellId_; }

private:
    Vector3d	   pos_;
    std::uint64_t  cellId_;
    double		   lum_;
    double         radius_;
};

struct SourceParameters
{
    std::uint64_t num_photons_;
    std::uint32_t PrimaryDirectionsLevel_;
    bool useMonteCarlo_;
    bool useHEALPixGrid_;
};

// all sources of photons
class Sources
{
public:
    Sources(SourceParameters parameters, std::vector<PointSource> pointSources, std::vector<SphereSource> sphereSources)
        : parameters_{ parameters}
        , pointSources_{ std::move(pointSources) }
        , sphereSources_{ std::move(sphereSources) }
        , totlum_{ 0. }
        , primaryDir_{ parameters_.useMonteCarlo_ ? 1 : parameters_.PrimaryDirectionsLevel_, parameters_.useHEALPixGrid_ }
        , sphereDir_{ 1 }
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

        for (const auto& sphereSource : sphereSources_)
        {
            totlum_ += sphereSource.luminosity();
        }

        double const firstLuminosity = pointSources_.empty()
            ? sphereSources_.at(0).luminosity()
            : pointSources_.at(0).luminosity();

        photonsNumber_ = parameters_.useMonteCarlo_
            ? (std::uint64_t) (parameters_.num_photons_ * firstLuminosity / totlum_)
            : primaryDir_.number();
    }

    ~Sources() = default;

    Sources(Sources const &) = delete;
    Sources& operator=(Sources const&) = delete;

    Photon emitPhoton(IGridCRef grid, IRandomGenerator* ran);
    void directPhotons(IGridCRef grid, std::vector<Observer>* observers);
    void writeObserversOpticalDepths(IGridCRef grid, std::vector<Observer>* observers);

    std::uint64_t num_photons() const
    {   return parameters_.num_photons_; }

private:
    SourceParameters parameters_;
    std::vector<PointSource> const pointSources_;
    std::vector<SphereSource> const sphereSources_;
    double	 totlum_;

    Directions primaryDir_;
    Directions sphereDir_;
    std::uint32_t currentSource_;
    std::uint64_t photonId_;
    std::uint64_t photonsNumber_;

    // dgem for sphere source
    std::uint64_t pointId_{0};
    Vector3d pointPosition_;
    std::uint64_t pointCellId_;
};

#endif
