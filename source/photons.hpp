#ifndef _PHOTONS_HPP_
#define _PHOTONS_HPP_

#include "model.hpp"

void Dustmat( double &p1, double &p2, double &p3, double &p4,
			double cost, double cost2, double pl, double pc, 
			double sc, double hgg, double g2 );

// photon. It's properties and methods for the work with it
class PHOTON
{
	public:
		// random direction photon generation
		PHOTON( POSITION const &pos, double weight, int nscat );
		// known direction photon generation
		PHOTON(POSITION const &pos, DIRECTION const &dir, double weight, int nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
		// photon scattering
		double Scatt( MODEL const &m, DIRECTION const & dir ); 
		void Scatt( MODEL const &m, DIRECTIONS const &dirs, GRID const &grid, DIRECTION const &obs, PICTURES *pict );
		void Stokes( MODEL const &m, DIRECTION const &dir, double calpha, bool fDir );
		void Move(double t)
		{
			pos_ = POSITION(pos_.x()+t*dir_.nx(), pos_.y()+t*dir_.ny(), pos_.z()+t*dir_.nz());
		}
		POSITION & pos() 
		{	return pos_;	}
		DIRECTION & dir() 
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
		POSITION	pos_;	// outpoint 
		DIRECTION	dir_;	// vector of the direction
		uint32_t 	nscat_;				// число рассеяний
		double weight_;				// статистический вес фотона
		double fi_, fq_, fu_, fv_;			// Stokes fluxes
};
// one source of photons
class SOURCE
{
	public:
		SOURCE() : pos_(POSITION(0.0, 0.0, 0.0)), lum_(0.0) {};
		SOURCE( POSITION const &pos, double lum ) : pos_(pos), lum_(lum) {};
		SOURCE( double x, double y, double z, double lum ) : pos_(POSITION(x,y,z)), lum_(lum) {};
		POSITION const & pos( void ) const
		{	return pos_; }
		double lum( void ) const
		{	return lum_; }
	private:
		POSITION	pos_;
		double		lum_; 
};
// all sources of photons
class SOURCES
{
	public:
		SOURCES(  ) : num_(0), totlum_(0), sources_(nullptr) {};
		~SOURCES()
		{
			if (num_ != 0) delete[] sources_;
		}

		void Init(uint32_t nstars, double *x, double *y, double *z, double *l)
		{
			num_ = nstars;
			sources_ = new SOURCE[num_];
			for (size_t cnt=0; cnt!=num_; ++cnt)
			{
				sources_[cnt] = SOURCE(x[cnt], y[cnt], z[cnt], l[cnt]);
				totlum_ += sources_[cnt].lum();
			}
		}
		uint32_t num( void ) const
		{	return num_;	}
		double totlum( void ) const
		{	return totlum_;	}
		SOURCE & operator []( int i ) 
		{	return sources_[i];	}
		SOURCE const & operator []( int i ) const
		{	return sources_[i];	}
	private:
		uint32_t 	num_;
		double	totlum_; 
		SOURCE	*sources_;

		SOURCES ( SOURCES const &);
		SOURCES & operator =( SOURCES const &);
};

// a structure for Stokes parameters in similar scatterings
class SCATHOLDER
{
	public:
		SCATHOLDER( bool fHold=0, double hgfac=0.0, DIRECTION const &dir=DIRECTION(0,0,0), double fi=0.0, double fq=0.0, double fu=0.0, double fv=0.0) :
			fHold_(fHold), hgfac_(hgfac), dir_(dir), fi_(fi), fq_(fq), fu_(fu), fv_(fv) {};
		bool fHold( void ) const
		{	return fHold_;	}
		double hgfac( void ) const
		{	return hgfac_;	}
		DIRECTION & dir( void ) 
		{	return dir_;	}
		double fi( void ) const
		{	return fi_;	}
		double fq( void ) const
		{	return fq_;	}
		double fu( void ) const
		{	return fu_;	}
		double fv( void ) const
		{	return fv_;	}	
	private:
		bool fHold_;
		double hgfac_;
		DIRECTION dir_;
		double fi_, fq_, fu_, fv_;
};
#endif
