#ifndef PHOTONS_HPP_
#define PHOTONS_HPP_

#include "model.hpp"

void Dustmat( double &p1, double &p2, double &p3, double &p4,
			double cost, double cost2, double pl, double pc, 
			double sc, double hgg, double g2 );

// photon. It's properties and methods for the work with it
class Photon
{
	public:
		// random direction photon generation
		Photon( Position const &pos, double weight, int nscat );
		// known direction photon generation
		Photon(Position const &pos, Direction const &dir, double weight, int nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
		// photon scattering
		double Scatt( Model const &m, Direction const & dir ); 
		void Scatt( Model const &m, Directions const &dirs, Grid const &grid, Direction const &obs, Pictures *pict );
		void Stokes( Model const &m, Direction const &dir, double calpha, bool fDir );
		void Move(double t)
		{
			pos_ = Position(pos_.x()+t*dir_.nx(), pos_.y()+t*dir_.ny(), pos_.z()+t*dir_.nz());
		}
		Position & pos() 
		{	return pos_;	}
		Direction & dir() 
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
		Position	pos_;	// outpoint 
		Direction	dir_;	// vector of the direction
		uint32_t 	nscat_;				// число рассеяний
		double weight_;				// статистический вес фотона
		double fi_, fq_, fu_, fv_;			// Stokes fluxes
};
// one source of photons
class Source
{
	public:
		Source() : pos_(Position(0.0, 0.0, 0.0)), lum_(0.0) {};
		Source( Position const &pos, double lum ) : pos_(pos), lum_(lum) {};
		Source( double x, double y, double z, double lum ) : pos_(Position(x,y,z)), lum_(lum) {};
		Position const & pos( void ) const
		{	return pos_; }
		double lum( void ) const
		{	return lum_; }
	private:
		Position	pos_;
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

// a structure for Stokes parameters in similar scatterings
class Scatholder
{
	public:
		Scatholder( bool fHold=0, double hgfac=0.0, Direction const &dir=Direction(0,0,0), double fi=0.0, double fq=0.0, double fu=0.0, double fv=0.0) :
			fHold_(fHold), hgfac_(hgfac), dir_(dir), fi_(fi), fq_(fq), fu_(fu), fv_(fv) {};
		bool fHold( void ) const
		{	return fHold_;	}
		double hgfac( void ) const
		{	return hgfac_;	}
		Direction & dir( void ) 
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
		Direction dir_;
		double fi_, fq_, fu_, fv_;
};
#endif
