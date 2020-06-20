#ifndef PHOTONS_HPP_
#define PHOTONS_HPP_

#include "model.hpp"
#include "Vector3d.hpp"
#include "Direction3d.hpp"

// photon. It's properties and methods for the work with it
class Photon
{
    public:
        // random direction photon generation
        Photon(Vector3d const &position, double weight, uint32_t nscat );
        // known direction photon generation
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

// all sources of photons
class Sources
{
    public:
        Sources(uint32_t nstars, double *x, double *y, double *z, double *l)
        : number_{ nstars }
        , totlum_{ 0. }
        {
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

        uint32_t num() const
        {	return number_;	}
        double totlum() const
        {	return totlum_;	}
        PointSource& operator []( int i )
        {	return sources_[i];	}
        PointSource const& operator []( int i ) const
        {	return sources_[i];	}
    private:
        uint32_t number_;
        double	 totlum_;
        PointSource	*sources_;
};

#endif
