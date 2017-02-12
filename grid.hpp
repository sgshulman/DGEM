#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <math.h>
#include "model.hpp"

class POSITION
{
	public:
		POSITION(double x=0.0, double y=0.0, double z=0.0) : x_(x), y_(y), z_(z)	{ };
		// data 
		double &x()
		{	return x_;	};	
		double &y()
		{	return y_;	};
		double &z()
		{	return z_;	};
	private:
		double x_, y_, z_;	
};

class DIRECTION
{
	public:
		// default zero init
		DIRECTION() : nx_(0.0), ny_(0.0), nz_(0.0), phi_(0.0), theta_(0.0), 
						sint_(0.0), cost_(0.0), sinp_(0.0), cosp_(0.0)
		{}; 
		DIRECTION(double nx, double ny, double nz)
		{
			double r=sqrt(nx*nx+ny*ny+nz*nz);			
			nx_ = nx/r; 
			ny_ = ny/r; 
			nz_ = nz/r;
			
			phi_ = atan2(ny_, nx_);
			if(phi_ > 2*3.1415926 ) phi_=phi_-2*3.1415926;
			if(phi_ < 0.0)          phi_=phi_+2*3.1415926;
			 
			theta_ = acos(nz_);

			cosp_ = cos(phi_);
			sinp_ = sin(phi_);
			sint_ = sin(theta_);
			cost_ = cos(theta_);
		}
		//
		DIRECTION(double phi, double theta)
		{
			phi_ = phi;
			if(phi_ > 2*3.1415926 ) phi_=phi_-2*3.1415926;
			if(phi_ < 0.0)          phi_=phi_+2*3.1415926;
			theta_ = theta;

			nx_ = sin(theta_)*cos(phi_);
			ny_ = sin(theta_)*sin(phi_); 
			nz_ = cos(theta_);

			cosp_ = cos(phi_);
			sinp_ = sin(phi_);
			sint_ = sin(theta_);
			cost_ = cos(theta_);
		}
		double nx() const
		{	return nx_;	}
		double ny() const
		{	return ny_;	}
		double nz() const
		{	return nz_;	}
		double phi() const
		{	return phi_;	}
		double cost() const
		{	return cost_;	}
		double sint() const
		{	return sint_;	}
		double cosp() const
		{	return cosp_;	}
		double sinp() const
		{	return sinp_;	}
	private:
		double nx_, ny_, nz_;
		double phi_, theta_;
		double sint_,cost_,sinp_,cosp_;
};
	
// cartezian grid
class GRID
{
	public:
		GRID() : Nx_(0), Ny_(0), Nz_(0), rhokappa_(nullptr), obstau_(nullptr) {};
		~GRID()
		{
			if (rhokappa_ != nullptr) delete[] rhokappa_;
			if (obstau_ != nullptr) delete[] obstau_;
		}
		void Init(const MODEL &m, double R_i, double R_d, double rho_0, double h_0, double R_0, 
									double alpha, double beta, uint32_t Nx, uint32_t Ny, uint32_t Nz );
		double PhotonSMax( PHOTON &ph ) const;
      	double PhotonCWall( PHOTON &ph, double delta ) const;
		double TauFind( PHOTON ph, double delta=-0.001 ) const;
		int TauInt( PHOTON & ph, double tau, double tauold=0.0, double delta=-0.001 ) const;
		int TauInt2( PHOTON & ph, double delta=-0.001 ) const;
		void Peeloff( PHOTON ph, DIRECTION const &obs, MODEL const &m, PICTURES *pict, SCATHOLDER *holder ) const;
		void Peeloff( PHOTON ph, MODEL const &m, PICTURES *pict, SCATHOLDER *holder ) const;
	private:
		uint32_t Nx_, Ny_, Nz_;
		double *rhokappa_;
		double *obstau_;	
		double xmax_, ymax_, zmax_;
		double minrho_;
		GRID ( GRID const &);
		GRID & operator =( GRID const &);
};
#endif

