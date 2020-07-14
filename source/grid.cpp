#include <algorithm>
#include <cmath>
#include <iostream>

#include "grid.hpp"
#include "observers.hpp"
#include "model.hpp"
#include "Photon.hpp"

namespace
{
    inline double ds(double const pos, double const velocity, double const max)
    {
        if(velocity > 0.0)
        {
            return (max - pos) / velocity;
        } else if(velocity < 0.0) {
            return -(pos + max) / velocity;
        }
        return 2.0 * max;
    }
}

double FlaredDisk::density(double x, double y, double z) const
{
    double const r2 = x*x + y*y;
    double const r  = std::sqrt(r2);

    // Disk Geometry
    if(( r >= rInner_ ) && ( r <= rOuter_ ))
    {
        double const h = h0_ * std::pow(r/r0_, beta_);
        return rho0_ * std::pow(r0_ / r, alpha_) * std::exp(-0.5*z*z / (h*h)); // rho in g/cm^3
    } else {
        return 0.0;
    }
}


Grid::Grid(
        double const xmax,
        double const ymax,
        double const zmax,
        double const kappa,
        uint32_t const nx,
        uint32_t const ny,
        uint32_t const nz,
        FlaredDiskCPtr disk)
    : nx_{ nx }
    , ny_{ ny }
    , nz_{ nz }
    , xmax_{ xmax }
    , ymax_{ ymax }
    , zmax_{ zmax }
    , disk_{ std::move(disk) }
{
    rhokappa_ = new double[nx_ * ny_ * nz_];
    minrho_ = 1e+38;

    for (size_t cntx=0; cntx!=nx_; ++cntx)
    {
        double const x = (cntx*2.0+1) * xmax_/nx_ - xmax_;
        for (size_t cnty=0; cnty!=ny_; ++cnty)
        {
            double const y=(cnty*2.0+1) * ymax_/ny_ - ymax_;
            for (size_t cntz=0; cntz!=nz_; ++cntz)
            {
                double const z=(cntz*2.0+1) * zmax_/nz_ - zmax_;
                size_t const idx = cntx+cnty*nx+cntz*ny*nx;
                rhokappa_[idx] = disk_->density(x, y, z) * kappa * 1.5e13; // rho*kappa*R,
                if (minrho_ > rhokappa_[idx] && rhokappa_[idx] > 0)
                    minrho_ = rhokappa_[idx];
            }
        }
    }
}

// calculate smax -- maximum distance photon can travel *******
double Grid::maxDistance(Photon const& ph) const
{	
    double const dsx = ds(ph.pos().x(), ph.dir().x(), xmax_);
    double const dsy = ds(ph.pos().y(), ph.dir().y(), ymax_);
    double const dsz = ds(ph.pos().z(), ph.dir().z(), zmax_);

    return std::min(dsx, std::min(dsy, dsz));
}

// find distance to next x, y, and z cell walls.  
// note that dx is not the x-distance, but the actual distance along 
// the direction of travel to the next x-face, and likewise for dy and dz.
double Grid::cellDistance(Photon& ph, double delta) const
{
    double dx=200.0*xmax_, dy=200.0*ymax_, dz=200.0*zmax_, dcell=0.0;
    if (delta < 0.0) delta=0.0001*(2.*xmax_/nx_);
    double tmp = 2.0*xmax_/nx_;
    if(ph.dir().x() > 0.0)
    {
        dx = (( int( (ph.pos().x()+xmax_) / tmp ) + 1.0 )* tmp - xmax_ - ph.pos().x() )/ph.dir().x();
        if (dx < delta)
        {
            dx = tmp/ph.dir().x();
            ph.pos().x() = (int( (ph.pos().x()+xmax_) / tmp ) + 1.0 )*tmp - xmax_ + delta;
        }
    } else if (ph.dir().x() < 0.0) {
        dx = (( int( (ph.pos().x()+xmax_)/ tmp ) )* tmp - xmax_ - ph.pos().x() )/ph.dir().x();
        if (dx < delta)
        {
            dx =-tmp/ph.dir().x();
            ph.pos().x() = ( int( (ph.pos().x()+xmax_)/tmp) )* tmp - xmax_ - delta;
        }
    }

    tmp = 2.0*ymax_/ny_;
    if(ph.dir().y() > 0.0)
    {
        dy = (( int( (ph.pos().y()+ymax_)/tmp ) + 1.0 )* tmp - ymax_ - ph.pos().y() )/ph.dir().y();
        if (dy < delta)
        {
            dy = tmp/ph.dir().y();
            ph.pos().y() = (int( (ph.pos().y()+ymax_)/tmp) + 1.0 )* tmp - ymax_+ delta;
        }
    } else if (ph.dir().y() < 0.0) {
        dy = (( int( (ph.pos().y()+ymax_)/tmp ) )*tmp - ymax_ - ph.pos().y() )/ph.dir().y();
        if (dy < delta)
        {
            dy =-tmp/ph.dir().y();
            ph.pos().y() = ( int( (ph.pos().y()+xmax_)/tmp ) )* tmp - ymax_ - delta;
        }
    }

    tmp=2.0*zmax_/nz_;
    if(ph.dir().z() > 0.0)
    {
        dz = (( int( (ph.pos().z()+zmax_)/tmp ) + 1.0 )* tmp - zmax_ - ph.pos().z() )/ph.dir().z();
        if (dz < delta)
        {
            dz = tmp/ph.dir().z();
            ph.pos().z() = (int( (ph.pos().z()+zmax_)/tmp ) + 1.0 )* tmp - zmax_ + delta;
        }
    } else if (ph.dir().z() < 0.0) {
        dz = (( int( (ph.pos().z()+zmax_)/tmp ))* tmp- zmax_ - ph.pos().z() )/ph.dir().z();
        if (dz < delta)
        {
            dz =-tmp/ph.dir().z();
            ph.pos().z() = ( int( (ph.pos().z()+zmax_)/tmp ) )* tmp - zmax_ - delta;
        }
    }

    dcell = ( dx   < dy ) ? dx    : dy;
    dcell = (dcell < dz ) ? dcell : dz;
//	printf("d %g %g %g \n", dx, dy, dz);
    return dcell;
}

double Grid::findOpticalDepth( Photon ph, double delta ) const
{
    double taurun=0.0, d=0.0;

    if (delta < 0.0) delta=0.0001*(2.*xmax_/nx_);
    double const smax = maxDistance(ph);

    if(smax < delta)
    {
        return 0.0;
    }

    while (d < 0.999*smax)
    {
        double const dcell = cellDistance(ph, delta);
        auto const x = uint32_t( (ph.pos().x()+xmax_)*nx_ / 2.0 / xmax_ );
        auto const y = uint32_t( (ph.pos().y()+ymax_)*ny_ / 2.0 / ymax_ );
        auto const z = uint32_t( (ph.pos().z()+zmax_)*nz_ / 2.0 / zmax_ );

        taurun += dcell * rhokappa_[ x+y*nx_+z*ny_*nx_ ];
        ph.Move(dcell);
        d += dcell;
    }

    return taurun;
}

int Grid::movePhotonAtDepth( Photon & ph, double tau, double tauold, double delta ) const
{
    double taurun=tauold, d=0.0;

    if (delta < 0.0 )
    {
        delta=0.0001*(2.*xmax_/nx_);
    }

    double const smax = maxDistance(ph);

    // integrate through grid
    while (taurun < tau && d < (0.9999*smax))
    {
        double const dcell = cellDistance(ph, delta);
        auto const x = uint32_t( (ph.pos().x()+xmax_)*nx_ / 2.0 / xmax_ );
        auto const y = uint32_t( (ph.pos().y()+ymax_)*ny_ / 2.0 / ymax_ );
        auto const z = uint32_t( (ph.pos().z()+zmax_)*nz_ / 2.0 / zmax_ );
        double const taucell = dcell*rhokappa_[ x+y*nx_+z*ny_*nx_ ];

        if(taurun + taucell >= tau)
        {
            double const d1 = (tau-taurun) / rhokappa_[ x+y*nx_+z*ny_*nx_ ];
            d += d1;
            ph.Move(d1);
        } else {
            d += dcell;
            ph.Move(dcell);
        }

        taurun += taucell;
    }

    // calculate photon final position.  if it escapes envelope then
    // set tflag=1.  if photon doesn't escape leave tflag=0 and update
    // photon position.
    return d >= 0.999 * smax;
}

int Grid::movePhotonAtRandomDepth( Photon &ph, Random *ran, double delta) const
{
    double const tau = -std::log(ran->Get());
    return movePhotonAtDepth(ph, tau, 0.0, delta);
}

void Grid::peeloff(Photon ph, Observer& observer, DustCRef dust) const
{
    double const hgfac = ph.Scatt(dust, observer.direction(), nullptr);
    double const tau = findOpticalDepth(ph);

    if (tau == 0.0)
    {
        return;
    }

    ph.weight() *= hgfac * exp(-tau);
    // Bin the photon into the image according to its position and direction of travel.
    observer.bin(ph);
} 
