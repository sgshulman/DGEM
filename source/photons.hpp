#ifndef PHOTONS_HPP_
#define PHOTONS_HPP_

#include <limits>
#include "model.hpp"
#include "Vector3d.hpp"
#include "Direction3d.hpp"
#include "directions.hpp"

// photon. It's properties and methods for the work with it
class Photon
{
    public:
        Photon(Vector3d const &position, Direction3d const &dir, double weight, uint32_t nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
        // photon scattering
        double Scatt(DustCRef dust, Direction3d const & dir );
        void Scatt( Model const &m, Directions const &dirs, GridCRef grid, std::vector<Observer>& observers);
        void Stokes(DustCRef dust, Direction3d const &dir, double calpha, bool fDir );
        void Move(double t)
        {
            pos_ = pos_ + t * dir_.vector();
        }
        Vector3d& pos()
        {	return pos_;	}
        Direction3d& dir()
        {	return dir_;	}
        Vector3d const& pos() const
        {	return pos_;	}
        Direction3d const& dir() const
        {	return dir_;	}

        double& x()
        {	return pos_.x();}
        double& y()
        {	return pos_.y();}
        double& z()
        {	return pos_.z();}

        double fi() const
        {	return fi_;		}
        double fq() const
        {	return fq_;		}
        double fu() const
        {	return fu_;		}
        double fv() const
        {	return fv_;		}

        uint32_t& nscat()
        {	return nscat_;	}
        double& weight()
        {	return weight_;	}

        bool termination() const
        {   return nscat_ == std::numeric_limits<uint32_t>::max(); }

    private:
        Vector3d	pos_;	// outpoint
        Direction3d	dir_;	// vector of the direction
        uint32_t 	nscat_;
        double weight_;
        double fi_, fq_, fu_, fv_;
};

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
            std::cout << "nstars " << nstars << std::endl;
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
