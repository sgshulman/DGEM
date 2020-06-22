#ifndef SOURCES_HPP_
#define SOURCES_HPP_

#include <cstdint>
#include "Vector3d.hpp"
#include "photons.hpp"
#include "Predefines.hpp"


// one point source of photons
class PointSource
{
public:
    PointSource()
        : lum_{}
    {}

    PointSource(Vector3d const& pos, double lum)
        : pos_(pos)
        , lum_(lum)
    {};

    Vector3d const& pos() const
    {	return pos_; }

    double luminosity() const
    {	return lum_; }

private:
    Vector3d	pos_;
    double		lum_;
};

struct SourceParameters
{
    bool useMonteCarlo_;
    uint64_t num_photons_;
    uint32_t PrimaryDirectionsLevel_;
    int32_t seed_;
};

// all sources of photons
class Sources
{
public:
    Sources(SourceParameters parameters, uint32_t nstars, double *x, double *y, double *z, double *l)
        : parameters_{ parameters}
        , number_{ nstars }
        , totlum_{ 0. }
        , random_{ parameters_.seed_ }
        , primaryDir_{ parameters_.PrimaryDirectionsLevel_ }
        , currentSource_{ 0 }
        , photonId_{ 0 }
    {
        if (!parameters.useMonteCarlo_)
        {
            parameters_.num_photons_ = primaryDir_.number();
        }

        sources_ = new PointSource[number_];
        for (uint32_t cnt=0; cnt!=number_; ++cnt)
        {
            sources_[cnt] = PointSource(Vector3d{x[cnt], y[cnt], z[cnt]}, l[cnt]);
            totlum_ += sources_[cnt].luminosity();
        }
    }

    ~Sources()
    {
        delete[] sources_;
    }

    Sources(Sources const &) = delete;
    Sources& operator=(Sources const&) = delete;

    Photon emitPhoton();
    void directPhotons(GridCRef grid, std::vector<Observer>* observers);

    size_t num_photons() const
    {   return parameters_.num_photons_; }

private:
    SourceParameters parameters_;
    uint32_t const number_;
    double	 totlum_;
    PointSource	*sources_;

    Random random_;
    Directions primaryDir_;
    uint32_t currentSource_;
    uint64_t photonId_;
};

#endif
