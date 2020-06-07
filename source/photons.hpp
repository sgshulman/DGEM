#ifndef PHOTONS_HPP_
#define PHOTONS_HPP_

#include "model.hpp"
#include "Vector3d.hpp"
#include "Direction3d.hpp"

void Dustmat( double &p1, double &p2, double &p3, double &p4,
            double cost, double cost2, double pl, double pc,
            double sc, double hgg, double g2 );

// photon. It's properties and methods for the work with it
class Photon
{
    public:
        // random direction photon generation
        Photon(Vector3d const &position, double weight, int nscat );
        // known direction photon generation
        Photon(Vector3d const &position, Direction3d const &dir, double weight, int nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
        // photon scattering
        double Scatt(std::shared_ptr<Dust const> const& dust, Direction3d const & dir );
        void Scatt( Model const &m, Directions const &dirs, Grid const &grid, std::vector<Observer>& observers);
        void Stokes(std::shared_ptr<Dust const> const& dust, Direction3d const &dir, double calpha, bool fDir );
        void Move(double t)
        {
            pos_ = pos_ + t * dir_.vector();
        }
        Vector3d& pos()
        {	return pos_;	}
        Direction3d& dir()
        {	return dir_;	}
        double &x( void )
        {	return pos_.x();}
        double &y( void )
        {	return pos_.y();}
        double &z( void )
        {	return pos_.z();}
        double fi( void ) const
        {	return fi_;		}
        double fq( void ) const
        {	return fq_;		}
        double fu( void ) const
        {	return fu_;		}
        double fv( void ) const
        {	return fv_;		}
        uint32_t & nscat( void )
        {	return nscat_;	}
        double & weight( void )
        {	return weight_;	}
    private:
        Vector3d	pos_;	// outpoint
        Direction3d	dir_;	// vector of the direction
        uint32_t 	nscat_;				// число рассеяний
        double weight_;				// статистический вес фотона
        double fi_, fq_, fu_, fv_;			// Stokes fluxes
};
// one source of photons
class Source
{
    public:
        Source() : pos_(Vector3d()), lum_(0.0) {};
        Source( Vector3d const& pos, double lum ) : pos_(pos), lum_(lum) {};
        Source( double x, double y, double z, double lum ) : pos_(Vector3d(x,y,z)), lum_(lum) {};
        Vector3d const & pos( void ) const
        {	return pos_; }
        double lum( void ) const
        {	return lum_; }
    private:
        Vector3d	pos_;
        double		lum_;
};
// all sources of photons
class Sources
{
    public:
        Sources(  ) : num_(0), totlum_(0), sources_(nullptr) {};
        ~Sources()
        {
            if (num_ != 0) delete[] sources_;
        }

        void Init(uint32_t nstars, double *x, double *y, double *z, double *l)
        {
            num_ = nstars;
            sources_ = new Source[num_];
            for (size_t cnt=0; cnt!=num_; ++cnt)
            {
                sources_[cnt] = Source(x[cnt], y[cnt], z[cnt], l[cnt]);
                totlum_ += sources_[cnt].lum();
            }
        }
        uint32_t num( void ) const
        {	return num_;	}
        double totlum( void ) const
        {	return totlum_;	}
        Source & operator []( int i )
        {	return sources_[i];	}
        Source const & operator []( int i ) const
        {	return sources_[i];	}
    private:
        uint32_t 	num_;
        double	totlum_;
        Source	*sources_;

        Sources ( Sources const &);
        Sources & operator =( Sources const &);
};

#endif
