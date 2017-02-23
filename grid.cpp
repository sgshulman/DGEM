#include <math.h>
#include "grid.hpp"
#include "inoutput.hpp"
#include "model.hpp"
#include "photons.hpp"
#include <iostream>

double Density(double x, double y, double z, double R_i, double R_d, double rho_0, double h_0, double R_0, double alpha, double beta)
{
	// model flat disk
	double r, r2, h;

	r2=x*x+y*y;
	r=sqrt(r2);

	// Disk Geometry 
	if(( r >= R_i ) && ( r <= R_d )) 
	{
		h=h_0*pow(r/R_0, beta);
		return (rho_0*pow(R_0/r, alpha)*exp(-0.5*z*z/(h*h))); // rho in g/cm^3
	} else {
        return 0.0;
	}
}


void GRID::Init(const MODEL &m, double R_i, double R_d, double rho_0, double h_0, double R_0, 
								double alpha, double beta, uint32_t Nx, uint32_t Ny, uint32_t Nz )
{
	double x, y, z;			
	Nx_=Nx;
	Ny_=Ny;
	Nz_=Nz;
	rhokappa_ = new double[Nx_*Ny_*Nz_];
	minrho_ = 1e+38;
	xmax_ = m.xmax();
	ymax_ = m.ymax();
	zmax_ = m.zmax();
	for (size_t cntx=0; cntx!=Nx_; ++cntx)
	{
		x=(cntx*2.0+1)*xmax_/Nx_-xmax_;
		for (size_t cnty=0; cnty!=Ny_; ++cnty)
		{
			y=(cnty*2.0+1)*ymax_/Ny_-ymax_;
			for (size_t cntz=0; cntz!=Nz_; ++cntz)
			{
				z=(cntz*2.0+1)*zmax_/Nz_-zmax_;
				rhokappa_[cntx+cnty*Nx+cntz*Ny*Nx]=Density(x,y,z, R_i, R_d, rho_0, h_0, R_0, alpha, beta)*m.kappa()*1.5e13; // rho*kappa*R,
				if (minrho_ > rhokappa_[cntx+cnty*Nx+cntz*Ny*Nx] && rhokappa_[cntx+cnty*Nx+cntz*Ny*Nx] > 0) 
					minrho_ = rhokappa_[cntx+cnty*Nx+cntz*Ny*Nx];
			}
		}
	}
}
// calculate smax -- maximum distance photon can travel *******
double GRID::PhotonSMax( PHOTON &ph ) const
{	
	double dsx=0.0, dsy=0.0, dsz=0.0, smax =0.0;
	if(ph.dir().nx() > 0.0) {
		dsx= ( xmax_ - ph.pos().x() )/ph.dir().nx();
	} else if(ph.dir().nx() < 0.0) {
		dsx=-( ph.pos().x() + xmax_ )/ph.dir().nx();
	} else {
		dsx=200.0*xmax_;
	}
	if(ph.dir().ny() > 0.0) {
		dsy= ( ymax_ - ph.pos().y() )/ph.dir().ny();
	} else if(ph.dir().ny() < 0.0) {
		dsy=-( ph.pos().y() + ymax_ )/ph.dir().ny();
	} else {
		dsy=200.0*ymax_;
	}
	if(ph.dir().nz() > 0.0) {
		dsz= ( zmax_ - ph.pos().z() )/ph.dir().nz();
	} else if(ph.dir().nz() < 0.0) {
		dsz=-( ph.pos().z() + zmax_ )/ph.dir().nz();
	} else {
		dsz=200.0*zmax_;
	}
	smax = ( dsx < dsy) ? dsx  : dsy;
	smax = (smax < dsz) ? smax : dsz;
	return smax;
}
// find distance to next x, y, and z cell walls.  
// note that dx is not the x-distance, but the actual distance along 
// the direction of travel to the next x-face, and likewise for dy and dz.
double GRID::PhotonCWall(PHOTON &ph, double delta) const
{
	double dx=200.0*xmax_, dy=200.0*ymax_, dz=200.0*zmax_, dcell=0.0;
	if (delta < 0.0 ) delta=0.0001*(2.*xmax_/Nx_);
	double tmp = 2.0*xmax_/Nx_;
	if(ph.dir().nx() > 0.0) 
	{
		dx = (( int( (ph.pos().x()+xmax_) / tmp ) + 1.0 )* tmp - xmax_ - ph.pos().x() )/ph.dir().nx();
		if (dx < delta)
		{
			dx = tmp/ph.dir().nx();
			ph.x() = (int( (ph.pos().x()+xmax_) / tmp ) + 1.0 )*tmp - xmax_ + delta;
		}
	} else if (ph.dir().nx() < 0.0) {
		dx = (( int( (ph.pos().x()+xmax_)/ tmp ) )* tmp - xmax_ - ph.pos().x() )/ph.dir().nx();
		if (dx < delta)		
		{
			dx =-tmp/ph.dir().nx();
			ph.x() = ( int( (ph.pos().x()+xmax_)/tmp) )* tmp - xmax_ - delta;
		}
	} 
	
	tmp = 2.0*ymax_/Ny_;
	if(ph.dir().ny() > 0.0) 
	{
		dy = (( int( (ph.pos().y()+ymax_)/tmp ) + 1.0 )* tmp - ymax_ - ph.pos().y() )/ph.dir().ny();
		if (dy < delta)		
		{
			dy = tmp/ph.dir().ny();
			ph.y() = (int( (ph.pos().y()+ymax_)/tmp) + 1.0 )* tmp - ymax_+ delta;
		}
	} else if (ph.dir().ny() < 0.0) {
		dy = (( int( (ph.pos().y()+ymax_)/tmp ) )*tmp - ymax_ - ph.pos().y() )/ph.dir().ny();
		if (dy < delta)		
		{
			dy =-tmp/ph.dir().ny();
			ph.y() = ( int( (ph.pos().y()+xmax_)/tmp ) )* tmp - ymax_ - delta;
		}
	} 
	
	tmp=2.0*zmax_/Nz_;
	if(ph.dir().nz() > 0.0) 
	{
		dz = (( int( (ph.pos().z()+zmax_)/tmp ) + 1.0 )* tmp - zmax_ - ph.pos().z() )/ph.dir().nz();
		if (dz < delta)		
		{
			dz = tmp/ph.dir().nz();
			ph.z() = (int( (ph.pos().z()+zmax_)/tmp ) + 1.0 )* tmp - zmax_ + delta;
		}
	} else if (ph.dir().nz() < 0.0) {
		dz = (( int( (ph.pos().z()+zmax_)/tmp ))* tmp- zmax_ - ph.pos().z() )/ph.dir().nz();
		if (dz < delta)		
		{
			dz =-tmp/ph.dir().nz();
			ph.z() = ( int( (ph.pos().z()+zmax_)/tmp ) )* tmp - zmax_ - delta;
		}
	}

	dcell = ( dx   < dy ) ? dx    : dy;
	dcell = (dcell < dz ) ? dcell : dz;
//	printf("d %g %g %g \n", dx, dy, dz);
	return dcell;
}

double GRID::TauFind( PHOTON ph, double delta ) const
{
	double taurun=0.0, taucell=0.0, d=0.0, dcell=0.0;
	double smax = 0.0;
	uint32_t x=0, y=0, z=0;
	if (delta < 0.0 ) delta=0.0001*(2.*xmax_/Nx_);
	smax = PhotonSMax(ph);
	
	if(smax < delta) return 0.0;
	while (d < 0.999*smax)
	{
		dcell =	PhotonCWall(ph, delta);
		x = uint32_t( (ph.pos().x()+xmax_)*Nx_ / 2.0 / xmax_ );
		y = uint32_t( (ph.pos().y()+ymax_)*Ny_ / 2.0 / ymax_ );
		z = uint32_t( (ph.pos().z()+zmax_)*Nz_ / 2.0 / zmax_ );		   
		   
		taucell=dcell*rhokappa_[ x+y*Nx_+z*Ny_*Nx_ ]; 
		taurun+=taucell;
		ph.Move( dcell );
		d += dcell;
	}
	return taurun;  
}

int GRID::TauInt( PHOTON & ph, double tau, double tauold, double delta ) const
{
	double taurun=tauold, taucell=0.0, d=0.0;
	double smax = 0.0, dcell = 0.0, d1=0.0;
	uint32_t x=0, y=0, z=0;
	if (delta < 0.0 ) delta=0.0001*(2.*xmax_/Nx_);
	smax = PhotonSMax(ph);
	// integrate through grid
	while ( taurun < tau && d < (0.9999*smax))
	{
		dcell =	PhotonCWall(ph, delta);
		
		x = uint32_t( (ph.pos().x()+xmax_)*Nx_ / 2.0 / xmax_ );
		y = uint32_t( (ph.pos().y()+ymax_)*Ny_ / 2.0 / ymax_ );
		z = uint32_t( (ph.pos().z()+zmax_)*Nz_ / 2.0 / zmax_ );
		taucell=dcell*rhokappa_[ x+y*Nx_+z*Ny_*Nx_ ]; // xcell,ycell,zcell
		if( (taurun+taucell) >= tau)
		{
			d1=(tau-taurun)/rhokappa_[ x+y*Nx_+z*Ny_*Nx_ ]; 
		    d=d+d1;
		    taurun=taurun+taucell;
		    ph.Move( d1 );
		} else {
		    d=d+dcell;
		    taurun=taurun+taucell;
		    ph.Move( dcell );
		}
	}
	if((d>=(0.999*smax))) return 1;
	return 0;
}

int GRID::TauInt2( PHOTON &ph, double delta ) const
{
	double tau=-log(ran.Get());
	double taurun=0.0, taucell=0.0, d=0.0, d1=0.0;
	double smax = 0.0, dcell = 0.0;
	uint32_t x, y, z;
	if (delta < 0.0 ) delta=0.0001*(2.*xmax_/Nx_);
	smax = PhotonSMax(ph);
	if(smax < delta) return 1;
	// integrate through grid
	while ((taurun < tau)&&(d<(0.9999*smax)))
	{
		dcell =	PhotonCWall(ph, delta);
		x = uint32_t( (ph.pos().x()+xmax_)*Nx_ / 2.0 / xmax_ );
		y = uint32_t( (ph.pos().y()+ymax_)*Ny_ / 2.0 / ymax_ );
		z = uint32_t( (ph.pos().z()+zmax_)*Nz_ / 2.0 / zmax_ );
		taucell=dcell*rhokappa_[ x+y*Nx_+z*Ny_*Nx_ ];
		if((taurun+taucell)>=tau) 
		{
			d1=(tau-taurun)/rhokappa_[ x+y*Nx_+z*Ny_*Nx_ ];
			d=d+d1;
			taurun=taurun+taucell;
			ph.Move( d1 );
		} else {
			d=d+dcell;
		    taurun=taurun+taucell;
		    ph.Move( dcell );
		}
	}
	// calculate photon final position.  if it escapes envelope then
	// set tflag=1.  if photon doesn't escape leave tflag=0 and update 
	// photon position.
	if((d>=(0.999*smax))) return 1;
	return 0;
}

void GRID::Peeloff( PHOTON ph, DIRECTION const & obs, MODEL const &m, PICTURES *pict, SCATHOLDER *holder ) const
{
	double hgfac = ph.Scatt( m, obs ); 
	
	if (holder != nullptr ) *holder = SCATHOLDER( true, hgfac, ph.dir(), ph.fi(), ph.fq(), ph.fu(), ph.fv() );
	
	double tau2 = TauFind(ph);
	
	if(tau2 == 0.0) return;
	double phot=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fi();
	double photq=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fq();
	double photu=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fu();
	// Bin the photon into the image according to its position and 
	//direction of travel. 
	pict[0].Bin( ph, phot, photq, photu);
	if ( ph.nscat() <= m.nscat() )
	{
		pict[ ph.nscat() ].Bin( ph, phot, photq, photu);
	}
} 
void GRID::Peeloff( PHOTON ph, MODEL const &m, PICTURES *pict,  SCATHOLDER *holder ) const
{
	double hgfac = holder->hgfac(); 
	
	ph = PHOTON( ph.pos(), holder->dir(), ph.weight(), ph.nscat(), holder->fi(), holder->fq(), holder->fu(), holder->fv() );
	
	double tau2 = TauFind(ph);

	if(tau2 == 0.0) return;
	double phot=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fi();
	double photq=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fq();
	double photu=ph.weight()*hgfac*m.albedo()*exp(-tau2)*ph.fu();
	// Bin the photon into the image according to its position and 
	//direction of travel. 
	pict[0].Bin( ph, phot, photq, photu);
	if ( ph.nscat() <= m.nscat() )
	{
		pict[ ph.nscat() ].Bin( ph, phot, photq, photu);
	}
}
