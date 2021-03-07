#include <algorithm>
#include <cmath>

#include "CartesianGrid.hpp"
#include "IMatter.hpp"
#include "observers.hpp"
#include "Photon.hpp"
#include "Random.hpp"
#include "Units.hpp"

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


CartesianGrid::CartesianGrid(
        double const xmax,
        double const ymax,
        double const zmax,
        double const kappa,
        std::uint32_t const nx,
        std::uint32_t const ny,
        std::uint32_t const nz,
        IMatterCPtr matter)
    : nx_{ nx }
    , ny_{ ny }
    , nz_{ nz }
    , xmax_{ xmax }
    , ymax_{ ymax }
    , zmax_{ zmax }
    , xCellSize_ { 2. * xmax_ / nx_ }
    , yCellSize_ { 2. * ymax_ / ny_ }
    , zCellSize_ { 2. * zmax_ / nz_ }
    , xCellSizeInv_{ 1. / xCellSize_ }
    , yCellSizeInv_{ 1. / yCellSize_ }
    , zCellSizeInv_{ 1. / zCellSize_ }
    , kappa_{ kappa }
    , matter_{ std::move(matter) }
{
    rhokappa_ = new double[nx_ * ny_ * nz_];
    minrho_ = 1e+38;

    for (std::uint64_t cntx=0; cntx!=nx_; ++cntx)
    {
        double const x = (cntx*2.0+1) * xmax_/nx_ - xmax_;
        for (std::uint64_t cnty=0; cnty!=ny_; ++cnty)
        {
            double const y=(cnty*2.0+1) * ymax_/ny_ - ymax_;
            for (std::uint64_t cntz=0; cntz!=nz_; ++cntz)
            {
                double const z=(cntz*2.0+1) * zmax_/nz_ - zmax_;
                std::uint64_t const idx = cntx+cnty*nx+cntz*ny*nx;
                rhokappa_[idx] = matter_->density({x, y, z}) * kappa * AU_Cm; // rho*kappa*R,
                if (minrho_ > rhokappa_[idx] && rhokappa_[idx] > 0)
                    minrho_ = rhokappa_[idx];
            }
        }
    }
}


double CartesianGrid::findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const
{
    double tau{ 0.0 };
    Vector3d pos{ position };
    Vector3d dirNormalized{ direction.normalized() };

    double const step{0.001 * std::min({1., xmax_, ymax_, zmax_})};

    while (std::abs(pos.x()) <= xmax_ && std::abs(pos.y()) <= ymax_ && std::abs(pos.z()) <= zmax_)
    {
        tau += step * matter_->density(pos) * kappa_ * AU_Cm;
        pos = pos + step * dirNormalized;
    }

    return tau;
}


// calculate smax -- maximum distance photon can travel *******
double CartesianGrid::maxDistance(Photon const& ph) const
{	
    double const dsx = ds(ph.pos().x(), ph.dir().x(), xmax_);
    double const dsy = ds(ph.pos().y(), ph.dir().y(), ymax_);
    double const dsz = ds(ph.pos().z(), ph.dir().z(), zmax_);

    return std::min(dsx, std::min(dsy, dsz));
}

// find distance to next x, y, and z cell walls.  
// note that dx is not the x-distance, but the actual distance along 
// the direction of travel to the next x-face, and likewise for dy and dz.
std::pair<double, std::uint64_t> CartesianGrid::cellDistance(Photon& ph, double delta, Vector3d const& phDirInv) const
{
    double dx=200.0*xmax_, dy=200.0*ymax_, dz=200.0*zmax_;
    std::uint64_t newCellX=0, newCellY=0, newCellZ=0;

    auto const x = static_cast<std::uint32_t>( ph.cellId() & 0x00000000FFFFu);
    auto const y = static_cast<std::uint32_t>((ph.cellId() & 0x0000FFFF0000u) >> 16u);
    auto const z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF00000000u) >> 32u);

    if(ph.dir().x() > 0.0)
    {
        dx = ((x + 1.0) * xCellSize_ - xmax_ - ph.pos().x()) * phDirInv.x();
        if (dx < delta)
        {
            dx = xCellSize_ * phDirInv.x();
            ph.pos().x() = (x + 1.0) * xCellSize_ - xmax_ + delta;
            ph.cellId() += 1u;
        }
        newCellX = ph.cellId() + 1u;
    } else if (ph.dir().x() < 0.0) {
        dx = (x * xCellSize_ - xmax_ - ph.pos().x()) * phDirInv.x();
        if (dx < delta)
        {
            dx = -xCellSize_ * phDirInv.x();
            ph.pos().x() = x * xCellSize_ - xmax_ - delta;
            ph.cellId() -= 1u;
        }
        newCellX = ph.cellId() - 1u;
    }

    if(ph.dir().y() > 0.0)
    {
        dy = ((y + 1.0) * yCellSize_ - ymax_ - ph.pos().y()) * phDirInv.y();
        if (dy < delta)
        {
            dy = yCellSize_ * phDirInv.y();
            ph.pos().y() = (y + 1.0 ) * yCellSize_ - ymax_+ delta;
            ph.cellId() += 0x000000010000u;
        }
        newCellY = ph.cellId() + 0x000000010000u;
    } else if (ph.dir().y() < 0.0) {
        dy = (y * yCellSize_ - ymax_ - ph.pos().y()) * phDirInv.y();
        if (dy < delta)
        {
            dy = -yCellSize_ * phDirInv.y();
            ph.pos().y() = y * yCellSize_ - ymax_ - delta;
            ph.cellId() -= 0x000000010000u;
        }
        newCellY = ph.cellId() - 0x000000010000u;
    }

    if(ph.dir().z() > 0.0)
    {
        dz = ((z + 1.0) * zCellSize_ - zmax_ - ph.pos().z()) * phDirInv.z();
        if (dz < delta)
        {
            dz = zCellSize_  * phDirInv.z();
            ph.pos().z() = (z + 1.0) * zCellSize_ - zmax_ + delta;
            ph.cellId() += 0x000100000000u;
        }
        newCellZ = ph.cellId() + 0x000100000000u;
    } else if (ph.dir().z() < 0.0) {
        dz = (z * zCellSize_- zmax_ - ph.pos().z()) * phDirInv.z();
        if (dz < delta)
        {
            dz = -zCellSize_ * phDirInv.z();
            ph.pos().z() = z * zCellSize_ - zmax_ - delta;
            ph.cellId() -= 0x000100000000u;
        }
        newCellZ = ph.cellId() - 0x000100000000u;
    }

    std::pair<double, std::uint64_t> dcell = (dx < dy) ? std::make_pair(dx, newCellX) : std::make_pair(dy, newCellY);
    dcell = (dcell.first < dz ) ? dcell : std::make_pair(dz, newCellZ);

    return dcell;
}

double CartesianGrid::findOpticalDepth(Photon ph) const
{
    double taurun=0.0, d=0.0;

    double const delta=0.0001 * xCellSize_;
    double const smax = maxDistance(ph);

    if(smax < delta)
    {
        return 0.0;
    }

    Vector3d const phDirInv = ph.dir().vector().inverse();

    while (d < 0.999*smax)
    {
        std::pair<double, std::uint64_t> const dcell = cellDistance(ph, delta, phDirInv);

        auto const x = static_cast<std::uint32_t>( ph.cellId() & 0x00000000FFFFu);
        auto const y = static_cast<std::uint32_t>((ph.cellId() & 0x0000FFFF0000u) >> 16u);
        auto const z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF00000000u) >> 32u);

        taurun += dcell.first * rhokappa_[ x+y*nx_+z*ny_*nx_ ];
        ph.Move(dcell.first, dcell.second);
        d += dcell.first;
    }

    return taurun;
}

int CartesianGrid::movePhotonAtDepth(Photon & ph, double tau, double tauold) const
{
    double taurun=tauold, d=0.0;

    double const delta=0.0001 * xCellSize_;
    double const smax = maxDistance(ph);
    Vector3d const phDirInv = ph.dir().vector().inverse();

    // integrate through grid
    while (taurun < tau && d < (0.9999*smax))
    {
        std::pair<double, std::uint64_t> const dcell = cellDistance(ph, delta, phDirInv);

        auto const x = static_cast<std::uint32_t>( ph.cellId() & 0x00000000FFFFu);
        auto const y = static_cast<std::uint32_t>((ph.cellId() & 0x0000FFFF0000u) >> 16u);
        auto const z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF00000000u) >> 32u);

        double const taucell = dcell.first * rhokappa_[ x+y*nx_+z*ny_*nx_ ];

        if(taurun + taucell >= tau)
        {
            double const d1 = (tau-taurun) / rhokappa_[ x+y*nx_+z*ny_*nx_ ];
            d += d1;
            ph.Move(d1, ph.cellId());
        } else {
            d += dcell.first;
            ph.Move(dcell.first, dcell.second);
        }

        taurun += taucell;
    }

    // calculate photon final position.  if it escapes envelope then
    // set tflag=1.  if photon doesn't escape leave tflag=0 and update
    // photon position.
    return d >= 0.999 * smax;
}


int CartesianGrid::movePhotonAtRandomDepth(Photon &ph, Random *ran) const
{
    double const tau = -std::log(ran->Get());
    return movePhotonAtDepth(ph, tau, 0.0);
}


void CartesianGrid::peeloff(Photon ph, Observer& observer, IDustCRef dust) const
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


double CartesianGrid::computeMatterMass() const
{
    double density = 0.0;

    for (std::uint64_t cntx=0; cntx!=nx_; ++cntx)
    {
        double const x = (cntx*2.0+1) * xmax_/nx_ - xmax_;
        for (std::uint64_t cnty=0; cnty!=ny_; ++cnty)
        {
            double const y=(cnty*2.0+1) * ymax_/ny_ - ymax_;
            for (std::uint64_t cntz=0; cntz!=nz_; ++cntz)
            {
                double const z=(cntz*2.0+1) * zmax_/nz_ - zmax_;
                density += matter_->density({x, y, z});
            }
        }
    }

    return density * GPerCm3_MSunPerAU3 * 8. * xmax_ / nx_ * ymax_ / ny_ * zmax_ / nz_;
}


std::uint64_t CartesianGrid::cellId(const Vector3d& position) const
{
    auto const x = static_cast<std::uint32_t>((position.x()+xmax_)*xCellSizeInv_);
    auto const y = static_cast<std::uint32_t>((position.y()+ymax_)*yCellSizeInv_);
    auto const z = static_cast<std::uint32_t>((position.z()+zmax_)*zCellSizeInv_);
    return x + y*0x000000010000u + z*0x000100000000u;
}
